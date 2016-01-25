import numpy as np

from math import acos, pi

from TriSoupMesh import TriSoupMesh
from Utilities import *

# A mesh composed of halfedge elements.
# This class follows a "fast and loose" python design philosophy. Data is stored
# on halfedges/vertices/edges/etc as it's created.
#   - staticGeometry=True means that the structue and positions of this mesh
#     will not be changed after creation. This means that geometry-derived values
#     will be cached internally for performance improvements.
class HalfEdgeMesh(object):

    ### Construct a halfedge mesh from a TriSoupMesh
    def __init__(self, soupMesh, readPosition=True, checkMesh=False, staticGeometry=True):

        ### Members

        # Sets of all of the non-fake objects. Note that these are somewhat sneakily named,
        # the do not include the imaginary halfedges/faces used to close boundaries. Those
        # are only tracked in the 'full' sets below, as the user generally shouldn't mess with
        # them
        self.halfEdges = set()
        self.verts = set()
        self.faces = set()
        self.edges = set()

        # These versions include imaginary objects.
        # TODO these names are lame
        self.halfEdgesFull = set()
        self.facesFull = set()

        print('\nConstructing HalfEdge Mesh...')

        # TODO typecheck to ensure the input is a soupMesh?

        # NOTE assumes faces have proper winding, may fail in bad ways otherwise.
        # TODO Detect bad things (LOTS of ways this could be broken)
        # TODO Recover from bad things


        # There are 3 steps for this process
        #   - Iterate through vertices, create vertex objects
        #   - Iterate through faces, creating face, edge, and halfedge objects
        #     (and connecting where possible)
        #   - Iterate through edges, connecting halfedge.twin

        # Find which vertices are actually used in the mesh
        usedVertInds = set()
        for f in soupMesh.tris:
            usedVertInds.add(f[0])
            usedVertInds.add(f[1])
            usedVertInds.add(f[2])
        nUnused = len(set([i for i in range(len(soupMesh.verts))]) - usedVertInds)
        if nUnused > 0:
            print('  Note: ' + str(nUnused) + ' vertices in the original mesh were not used in any face and are being discarded')

        # Create vertex objects for only the used verts
        verts = []
        for (i, soupVert) in enumerate(soupMesh.verts):
            if i in usedVertInds:
                if readPosition:
                    v = Vertex(soupVert, staticGeometry=staticGeometry)
                else:
                    v = Vertex(staticGeometry=staticGeometry)
                verts.append(v)
                self.verts.add(v)
            else:
                verts.append(None)

        # Iterate over the faces, creating a new face and new edges/halfedges
        # for each. Fill out all properties except the twin & next references
        # for the halfedge, which will be handled below.
        edgeDict = {}           # The edge that connects two verts  [(ind1, ind2) ==> edge]
        edgeHalfEdgeDict = {}   # The two halfedges that border an edge [(ind1, ind2) ==> [halfedge list]]
        edgeSet = set()         # All edges that appear in the mesh, used for a sanity check
        for soupFace in soupMesh.tris:

            face = Face(staticGeometry=staticGeometry)
            self.faces.add(face)
            self.facesFull.add(face)

            theseHalfEdges = []     # The halfedges that make up this face

            # Iterate over the edges that make up the face
            for i in range(3):

                ind1 = soupFace[i]
                ind2 = soupFace[(i+1)%3]
                edgeKey = tuple(sorted((ind1,ind2)))

                # Sanity check that there are no duplicate edges in the input mesh
                if((ind1, ind2) in edgeSet):
                    raise ValueError('Mesh has duplicate edges or inconsistent winding, cannot represent as a half-edge mesh')
                else:
                    edgeSet.add((ind1, ind2))

                # Get an edge object, creating a new one if needed
                if(edgeKey in edgeDict):
                    edge = edgeDict[edgeKey]
                else:
                    edge = Edge(staticGeometry=staticGeometry)
                    edgeDict[edgeKey] = edge
                    self.edges.add(edge)
                    edgeHalfEdgeDict[edgeKey] = []

                # Create a new halfedge, which is always needed
                h = HalfEdge(staticGeometry=staticGeometry)
                self.halfEdges.add(h)
                self.halfEdgesFull.add(h)
                theseHalfEdges.append(h)

                # Set references to the halfedge in the other structures
                # This might be overwriting a previous value, but that's fine
                face.anyHalfEdge = h
                edge.anyHalfEdge = h
                verts[ind1].anyHalfEdge = h

                edgeHalfEdgeDict[edgeKey].append(h)

                # Set references to the other structures in the halfedge
                h.vertex = verts[ind2]
                h.edge = edge
                h.face = face


            # Connect the halfEdge.next reference for each of the halfedges we just created
            # in this face
            for i in range(3):
                theseHalfEdges[i].nextHE = theseHalfEdges[(i+1)%3]

        # Sanity check on the edges we say
        unpairedEdges = 0
        unpairedVerts = set()
        for (v1, v2) in edgeSet:
            if (v2, v1) not in edgeSet:
                unpairedEdges += 1
                unpairedVerts.add(v1)
                unpairedVerts.add(v2)
        print('  Input mesh has ' + str(unpairedEdges) + ' unpaired edges (which only appear in one direction)')
        print('  Input mesh has ' + str(len(unpairedVerts)) + ' unpaired verts (which touch some unpaired edge)')



        # Iterate through the edges to fill out the twin reference for each halfedge
        # This is where we use edgeHalfEdgeDict.
        for (edgeKey, halfEdgeList) in edgeHalfEdgeDict.iteritems():

            # print(edgeKey)
            # print(halfEdgeList)

            # Assuming the mesh is well-formed, this must be a list with two elements
            if(len(halfEdgeList) == 2):
                halfEdgeList[0].twinHE = halfEdgeList[1]
                halfEdgeList[1].twinHE = halfEdgeList[0]
            elif(len(halfEdgeList) > 2):
                raise ValueError('Mesh has more than two faces meeting at some edge')


        # Close boundaries by iterating around each hole and creating an imaginary face to cover the hole,
        # along with the associated halfedges. Note that this face will not be a triangle, in general.
        initialHalfEdges = self.halfEdges.copy()
        nHolesFilled = 0
        for initialHE in initialHalfEdges:

            # If this halfedge has no twin, then we have found a new boundary hole. Traverse the outside
            # and create a new faces/new halfedges
            # Note: strange things will happen if the multiples holes touch a single vertex.
            if initialHE.twinHE is None:
                nHolesFilled += 1

                fakeFace = Face(isReal=False, staticGeometry=staticGeometry)
                self.facesFull.add(fakeFace)

                # Traverse around the outside of the hole
                currRealHE = initialHE
                prevNewHE = None
                while True:

                    # Create a new fake halfedge
                    currNewHE = HalfEdge(isReal=False, staticGeometry=staticGeometry)
                    self.halfEdgesFull.add(currNewHE)
                    currNewHE.twinHE = currRealHE
                    currRealHE.twinHE = currNewHE
                    currNewHE.face = fakeFace
                    currNewHE.vertex = currRealHE.nextHE.nextHE.vertex
                    currNewHE.edge = currRealHE.edge
                    currNewHE.nextHE = prevNewHE

                    # Advance to the next border vertex along the loop
                    currRealHE = currRealHE.nextHE
                    while currRealHE != initialHE and currRealHE.twinHE != None:
                        currRealHE = currRealHE.twinHE.nextHE

                    prevNewHE = currNewHE

                    # Terminate when we have walked all the way around the loop
                    if currRealHE == initialHE:
                        break


                # Arbitrary point the fakeFace at the last created halfedge
                fakeFace.anyHalfEdge = currNewHE

                # Connect the nextHE ref for the first face edge, which was missed in the above loop
                initialHE.twinHE.nextHE = prevNewHE

        print('  Filled %d boundary holes in mesh using imaginary halfedges/faces'%(nHolesFilled))


        print("HalfEdge mesh construction completed")

        # Print out statistics about the mesh and check it
        self.printMeshStats(printImaginary=True)
        if checkMesh:
            self.checkMeshReferences()
        #self.checkDegenerateFaces() # a lot of meshes fail this...

    # Perform a basic refence validity check to catch blatant errors
    # Throws and error if it finds something broken about the datastructure
    def checkMeshReferences(self):

        # TODO why does this use AssertionError() instead of just assert statements?
        print('Testing mesh for obvious problems...')

        # Make sure the 'full' sets are a subset of their non-full counterparts
        diff = self.halfEdges - self.halfEdgesFull
        if(diff):
            raise AssertionError('ERROR: Mesh check failed. halfEdges is not a subset of halfEdgesFull')
        diff = self.faces - self.facesFull
        if(diff):
            raise AssertionError('ERROR: Mesh check failed. faces is not a subset of facesFull')

        # Accumulators for things that were referenced somewhere
        allRefHalfEdges = set()
        allRefEdges = set()
        allRefFaces = set()
        allRefVerts= set()

        ## Verify that every object in our sets is referenced by some halfedge, and vice versa
        # Accumulate sets of anything referenced anywhere and ensure no references are None
        for he in self.halfEdgesFull:

            if not he.nextHE:
                raise AssertionError('ERROR: Mesh check failed. he.nextHE is None')
            if not he.twinHE:
                raise AssertionError('ERROR: Mesh check failed. he.twinHE is None')
            if not he.edge:
                raise AssertionError('ERROR: Mesh check failed. he.edge is None')
            if not he.face:
                raise AssertionError('ERROR: Mesh check failed. he.face is None')
            if not he.vertex:
                raise AssertionError('ERROR: Mesh check failed. he.vertex is None')

            allRefHalfEdges.add(he.nextHE)
            allRefHalfEdges.add(he.twinHE)
            allRefEdges.add(he.edge)
            allRefFaces.add(he.face)
            allRefVerts.add(he.vertex)

            if he.twinHE.twinHE != he:
                raise AssertionError('ERROR: Mesh check failed. he.twin symmetry broken')

        for edge in self.edges:

            if not edge.anyHalfEdge:
                raise AssertionError('ERROR: Mesh check failed. edge.anyHalfEdge is None')
            allRefHalfEdges.add(edge.anyHalfEdge)

        for vert in self.verts:

            if not vert.anyHalfEdge:
                raise AssertionError('ERROR: Mesh check failed. vert.anyHalfEdge is None')
            allRefHalfEdges.add(vert.anyHalfEdge)

        for face in self.facesFull:

            if not face.anyHalfEdge:
                raise AssertionError('ERROR: Mesh check failed. face.anyHalfEdge is None')
            allRefHalfEdges.add(face.anyHalfEdge)

        # Check the resulting sets for equality
        if allRefHalfEdges != self.halfEdgesFull:
            raise AssertionError('ERROR: Mesh check failed. Referenced halfedges do not match halfedge set')
        if allRefEdges != self.edges:
            raise AssertionError('ERROR: Mesh check failed. Referenced edges do not match edges set')
        if allRefFaces != self.facesFull:
            raise AssertionError('ERROR: Mesh check failed. Referenced faces do not match faces set')
        if allRefVerts != self.verts:
            raise AssertionError('ERROR: Mesh check failed. Referenced verts do not match verts set')

        print('  ...test passed')


    def checkDegenerateFaces(self):
        """
        Checks if the mesh has any degenerate faces, which can mess up many algorithms.
        This is an exact-comparison check, so it won't catch vertices that differ by epsilon.
        """
        print("Checking mesh for degenerate faces...")

        for face in self.faces:

            seenPos = set()
            vList = []
            for v in face.adjVerts():
                pos = tuple(v.pos.tolist()) # need it as a hashable type
                if pos in seenPos:
                    raise ValueError("ERROR: Degenerate mesh face has repeated vertices at position: " + str(pos))
                else:
                    seenPos.add(pos)
                vList.append(v.pos)

            # Check for triangular faces with colinear vertices (don't catch other such errors for now)
            if(len(vList) == 3):
                v1 = vList[1] - vList[0]
                v2 = vList[2]-vList[0]
                area = norm(cross(v1, v2))
                if area < 0.0000000001*max((norm(v1),norm(v2))):
                    raise ValueError("ERROR: Degenerate mesh face has triangle composed of 3 colinear points: \
                        " + str(vList))


        print("  ...test passed")

    # Print out some summary statistics about the mesh
    def printMeshStats(self, printImaginary=False):

        if printImaginary:
            print('=== HalfEdge mesh statistics:')
            print('    Halfedges = %d  (+ %d imaginary)'%(len(self.halfEdges), (len(self.halfEdgesFull) - len(self.halfEdges))))
            print('    Edges = %d'%(len(self.edges)))
            print('    Faces = %d  (+ %d imaginary)'%(len(self.faces), (len(self.facesFull) - len(self.faces))))
            print('    Verts = %d'%(len(self.verts)))
        else:
            print('=== HalfEdge mesh statistics:')
            print('    Halfedges = %d'%(len(self.halfEdges)))
            print('    Edges = %d'%(len(self.edges)))
            print('    Faces = %d'%(len(self.faces)))
            print('    Verts = %d'%(len(self.verts)))


        maxDegree = max([v.degree() for v in self.verts])
        minDegree = min([v.degree() for v in self.verts])
        print('    - Max vertex degree = ' + str(maxDegree))
        print('    - Min vertex degree = ' + str(minDegree))

        nBoundaryVerts = sum([v.isBoundary() for v in self.verts])
        print('    - n boundary verts = ' + str(nBoundaryVerts))


    # If this mesh has boundaries, close them (make the imaginary faces/halfedges real).
    # The resulting mesh will be manifold.
    # Note: This naively triangulates non-triangular boundary faces, and thus can create
    # low quality meshes
    # TODO implement
    def fillBoundaries(self):
        raise NotImplementedError('fillBoundaries is not yet implemented')


    def enumerateVertices(self, subset=None):
        """
        Return a dictionary which assigns a 0-indexed integer to each vertex
        in the mesh. If 'subset' is given (should be a set), only the vertices
        in subset are indexed.
        """
        if subset is None:
            subset = self.verts

        enum = dict()
        ind = 0
        for vert in subset:
            if vert not in self.verts:
                raise ValueError("ERROR: enumerateVertices(subset) was called with a vertex in subset which is not in the mesh.")

            enum[vert] = ind
            ind += 1

        return enum


    def assignReferenceDirections(self):
        '''
        For each vertex in the mesh, arbitrarily selects one outgoing halfedge
        as a reference ('refEdge').
        '''
        for vert in self.verts:
            vert.refEdge = vert.anyHalfEdge


    def applyVertexValue(self, value, attributeName):
        """
        Given a dictionary of {vertex => value}, stores that value on each vertex
        with attributeName
        """

        # Throw an error if there isn't a value for every vertex
        if not set(value.keys()) == self.verts:
            raise ValueError("ERROR: Attempted to apply vertex values from a map whos domain is not the vertex set")

        for v in self.verts:
            setattr(v, attributeName, value[v])

    # Returns a brand new TriSoupMesh corresponding to this mesh
    # 'retainVertAttr' is a list of vertex-valued attributes to carry in to th trisoupmesh
    # TODO do face attributes (and maybe edge?)
    # TODO Maybe implement a 'view' version of this, so that we can modify the HalfEdge mesh
    # without completely recreating a new TriSoup mesh.
    def toTriSoupmesh(self,retainVertAttr=[]):

        # Create a dictionary for the vertex attributes we will retain
        vertAttr = dict()
        for attr in retainVertAttr:
            vertAttr[attr] = []

        # Iterate over the vertices, numbering them and building an array
        vertArr = []
        vertInd = {}
        for (ind, v) in enumerate(self.verts):
            vertArr.append(v.pos)
            vertInd[v] = ind

            # Add any vertex attributes to the list
            for attr in retainVertAttr:
                vertAttr[attr].append(getattr(v, attr))

        # Iterate over the faces, building a list of the verts for each
        faces = []
        for face in self.faces:

            # Get the three edges which make up this face
            he1 = face.anyHalfEdge
            he2 = he1.nextHE
            he3 = he2.nextHE

            # Get the three vertices that make up the face
            v1 = vertInd[he1.vertex]
            v2 = vertInd[he2.vertex]
            v3 = vertInd[he3.vertex]

            faceInd = [v1, v2, v3]
            faces.append(faceInd)


        soupMesh = TriSoupMesh(vertArr, faces, vertAttr=vertAttr)

        return soupMesh


    # Computes face normals and area-weighted vertex normals.
    # Free bonus: computes face areas
    # def computeNormals(self):
    #
    #     print("Computing normals...")
    #
    #     # For each face, compute its normal and area
    #     for face in self.faces:
    #
    #         v = list(face.adjVerts())
    #
    #         faceNormal = np.cross(v[1].pos - v[0].pos, v[2].pos - v[0].pos)
    #         face.area = 0.5 * norm(faceNormal)
    #         face.normal = normalize(faceNormal)
    #
    #     # For each vertex, let its normal be the area-weighted average of
    #     # the adjacent faces
    #     for vert in self.verts:
    #         normalSum = np.array([0.0,0.0,0.0])
    #
    #         for face in vert.adjFaces():
    #             normalSum += face.normal * face.area
    #
    #         vert.normal = normalize(normalSum)
    #
    #
    #     print("Done computing normals")


class HalfEdge(object):

    ### Construct a halfedge, possibly not real
    def __init__(self, isReal=True, staticGeometry=False):
        self.isReal = isReal  # Is this a real halfedge, or an imaginary one we created to close a boundary?

        ### Members
        self.twinHE = None
        self.nextHE = None   # TODO this is a stupid name but 'next' is reserved
        self.vertex = None
        self.edge = None
        self.face = None

        self._cache = dict()
        self.staticGeometry = staticGeometry


    # Return a boolean indicating whether this is on the boundary of the mesh
    def isBoundary(self):
        return not self.twinHE.isReal

    @property
    def vec(self):
        """The vector represented by this halfedge"""
        if 'vec' in self._cache: return self._cache['vec']
        v = self.vertex.pos - self.twinHE.vertex.pos
        if self.staticGeometry: self._cache['vec'] = v
        return v


    @property
    def cotan(self):
        """
        Return the cotangent of the opposite angle, or 0 if this is an imaginary
        halfedge
        """
        # Validate that this is on a triangle
        if self.nextHE.nextHE.nextHE is not self:
            raise ValueError("ERROR: halfedge.cotan() is only well-defined on a triangle")

        if self.isReal:
            if 'cotan' in self._cache: return self._cache['cotan']

            # TODO implement me!
            val = 0.0

            if self.staticGeometry: self._cache['cotan'] = val
            return val

        else:
            return 0.0


    @property
    def oppAngle(self):
        """
        Compute the angle between the vectors opposite this edge in a triangle.
        Only well-defined on a triangle mesh.
        """
        if 'oppAngle' in self._cache: return self._cache['oppAngle']

        # Validate that this is on a triangle
        if self.nextHE.nextHE.nextHE is not self:
            print("isReal = " + str(self.isReal))
            raise ValueError("ERROR: halfedge.oppAngle() is only well-defined on a triangle")

        A = -normalize(self.nextHE.vec)
        B = normalize(self.nextHE.nextHE.vec)
        oppAngle = acos(np.dot(A,B))

        # NOTE: No triangle can have an angle greater than pi, so we don't need
        # to worry about correct for the special case of an angle greater than
        # pi as we do elsewhere

        if self.staticGeometry: self._cache['oppAngle'] = oppAngle
        return oppAngle

    @property
    def rescaledOppAngle(self):
        """
        The opposige angle (see oppAngle), rescaled so the vertex would
        have no angle defect.
        """
        if 'rescaledOppAngle' in self._cache: return self._cache['rescaledOppAngle']

        centerVert = self.twinHE.vertex
        angSum = 2*pi + centerVert.angleDefect
        scaledAngle = 2*pi * self.oppAngle / angSum

        if self.staticGeometry: self._cache['rescaledOppAngle'] = scaledAngle
        return scaledAngle

    @property
    def edgeAngle(self):
        """
        The angle coordinate (in radians) for this edge against the arbitrary reference
        direction defined in the tangent space of the source vertex
        """
        # NOTE: This method works by summing around the outgoing halfedges for the
        # vertex, which means that calling this as a property could be significantly
        # slower than pre-computing all of the values.

        if 'edgeAngle' in self._cache: return self._cache['edgeAngle']

        # Start at the reference edge and sum up angles until we reach this edge
        sourceVert = self.twinHE.vertex
        angSum = 0.0
        currOutEdge = sourceVert.refEdge
        while currOutEdge is not self:
            # TODO can optimize this loop
            nextEdge = currOutEdge.twinHE.nextHE
            A = normalize(currOutEdge.vec)
            B = normalize(nextEdge.vec)
            angSum += acos(np.dot(A,B))

            # If the vectors make an angle of more than pi, we would be choosing
            # the wrong inverse of cosine. Need to detect and correct for this case
            if np.dot(cross(B,A), sourceVert.normal) < 0:
                angSum += pi

            currOutEdge = nextEdge
        # print("    Raw angsum = " + str(angSum))
        # Since the loop above winds the wrong way, we need to subract to
        # to get the actual result measured CCW
        totalAngSum = 2*pi + sourceVert.angleDefect
        # print("    Total angsum = " + str(totalAngSum))
        angle = totalAngSum - angSum

        if self.staticGeometry: self._cache['edgeAngle'] = angle
        return angle


    @property
    def edgeAngleRescaled(self):
        """
        The angle coordinate (in radians) for this edge against the arbitrary reference
        direction defined in the tangent space of the source vertex, rescaled for the
        angle defect so that the angle sum is 2pi.
        """

        if 'edgeAngleRescaled' in self._cache:  return self._cache['edgeAngleRescaled']

        sourceVert = self.twinHE.vertex
        vertAngleSum = sourceVert.angleDefect + 2*pi
        scaledAngle = 2*pi * self.edgeAngle / vertAngleSum

        if self.staticGeometry: self._cache['edgeAngleRescaled'] = scaledAngle
        # print("  Rescaled edge angle is: " + str(scaledAngle))
        return scaledAngle

    @property
    def transportAngle(self):
        """
        The rotation (in radians) for parallel transport along this edge
        """
        if 'transportAngle' in self._cache: return self._cache['transportAngle']

        angle = regAngle(self.twinHE.edgeAngleRescaled - self.edgeAngleRescaled + pi)

        # print("\nTransport angle")
        # print("  vSource has ref " + str(self.twinHE.vertex.refEdge.vec))
        # print("  vDest   has ref " + str(self.vertex.refEdge.vec))
        # print("  vDest   has angle " + str(self.twinHE.edgeAngleRescaled))
        # print("  vSource has angle " + str(self.edgeAngleRescaled))
        # print("  transport angle is " + str(angle))

        if self.staticGeometry: self._cache['transportAngle'] = angle
        # print("Transport angle is: " + str(angle))
        return angle


class Vertex(object):

    ### Construct a vertex, possibly with a known position
    def __init__(self, pos=None, staticGeometry=False):

        if pos is not None:
            self._pos = pos
            if staticGeometry:
                self._pos.flags.writeable = False

        self.anyHalfEdge = None      # Any halfedge exiting this vertex

        self._cache = dict()
        self.staticGeometry = staticGeometry

    @property
    def pos(self):
        return self._pos

    @pos.setter
    def pos(self, value):
        if self.staticGeometry:
            raise ValueError("ERROR: Cannot write to vertex position with staticGeometry=True. To allow dynamic geometry, set staticGeometry=False when creating vertex (or in the parent mesh constructor)")
        self._pos = value


    # Iterate over the faces adjacent to this vertex (skips imaginary faces by default)
    def adjFaces(self, skipImaginary=True):

        # Iterate through the adjacent faces
        first = self.anyHalfEdge
        curr = self.anyHalfEdge
        while True:

            # Yield only real faces
            if curr.isReal or not skipImaginary:
                yield curr.face

            curr = curr.twinHE.nextHE
            if(curr == first):
                break

    # Iterate over the edges adjacent to this vertex
    def adjEdges(self):

        # Iterate through the adjacent edges
        first = self.anyHalfEdge
        curr = self.anyHalfEdge
        while True:

            yield curr.edge

            curr = curr.twinHE.nextHE
            if(curr == first):
                break
    # Iterate over the halfedges adjacent to this vertex
    def adjHalfEdges(self):

        # Iterate through the adjacent edges
        first = self.anyHalfEdge
        curr = self.anyHalfEdge
        while True:

            yield curr

            curr = curr.twinHE.nextHE
            if(curr == first):
                break

    # Iterate over the verts adjacent to this vertex
    def adjVerts(self):

        # Iterate through the adjacent edges
        first = self.anyHalfEdge
        curr = self.anyHalfEdge
        while True:

            yield curr.vertex

            curr = curr.twinHE.nextHE
            if(curr == first):
                break

    def adjEdgeVerts(self):
        """
        Iterate through the neighbors of this vertex, yielding a (edge,vert) tuple
        """
        first = self.anyHalfEdge
        curr = self.anyHalfEdge
        while True:

            yield (curr.edge, curr.vertex)

            curr = curr.twinHE.nextHE
            if(curr == first):
                break


    # Return a boolean indicating whether this is on the boundary of the mesh
    def isBoundary(self):

        # Traverse the halfedges adjacent to this, a loop of non-boundary halfedges
        # indicates that this vert is internal
        first = self.anyHalfEdge
        curr = self.anyHalfEdge
        while True:
            if curr.isBoundary():
                return True

            curr = curr.twinHE.nextHE
            if(curr == first):
                break

        return False

    # Returns the number edges/faces neighboring this vertex
    # TODO make this and similars use properties
    def degree(self):
        if 'degree' in self._cache: return self._cache['degree']

        d = sum(1 for e in self.adjEdges())

        if self.staticGeometry: self._cache['degree'] = d
        return d


    @property
    def normal(self):
        """The area-weighted normal vector for this face"""
        if 'normal' in self._cache: return self._cache['normal']

        normalSum = np.array([0.0,0.0,0.0])
        for face in self.adjFaces():
            normalSum += face.normal * face.area
        n = normalize(normalSum)

        if self.staticGeometry: self._cache['normal'] = n
        return n


    @property
    def angleDefect(self):
        """
        Compute the angle defect, the difference between the sum of the wedge-angles
        for all triangles adjacent to this and 2pi (return S - 2pi).
        For a boundary vertex, this function defines the curvature defect to be
        zero.
        """
        if 'angleDefect' in self._cache: return self._cache['angleDefect']

        if(self.isBoundary()):
            defect = 0

        else:

            angSum = 0.0
            vList = [normalize(h.vec) for h in self.adjHalfEdges()]
            # print("    " + str(vList))
            for (a,b) in circlePairs(vList):
                oppAngle = acos(np.dot(a,b))

                # If the vectors make an angle of more than pi, we would be choosing
                # the wrong inverse of cosine. Need to detect and correct for this case
                if np.dot(cross(b,a), self.normal) < 0:
                    oppAngle += pi

                # print("      + " + str(oppAngle))
                angSum += oppAngle

            defect = angSum - 2*pi

        if self.staticGeometry: self._cache['angleDefect'] = defect
        return defect

    @property
    def refDirR3(self):
        """
        Computes a reference direction ('refDirR3') in R3, guaranteed to
        lie in the tangent plane (as defined by the normal). This direction is
        the projection of the reference edge on to the tangent plane.

        Should only be called after a reference edge has been selected.
        (mesh.assignReferenceDirections())
        """
        if 'refDirR3' in self._cache: return self._cache['refDirR3']

        r = self.refEdge.vec
        rNorm = np.dot(self.normal, r)*self.normal
        rInPlane = normalize(r - rNorm)
        refDirR3 = rInPlane

        if self.staticGeometry: self._cache['refDirR3'] = refDirR3
        return refDirR3

class Face(object):


    ### Construct a face, possibly not real
    def __init__(self, isReal=True, staticGeometry=False):

        ### Members
        self.anyHalfEdge = None      # Any halfedge bordering this face
        self.isReal = isReal         # Is this an actual face of the mesh, or an artificial face we
                                     # created to close boundaries?

        self._cache = dict()
        self.staticGeometry = staticGeometry

    # Return a boolean indicating whether this is on the boundary of the mesh
    def isBoundary(self):

        first = self.anyHalfEdge
        curr = self.anyHalfEdge
        while True:
            if curr.isBoundary():
                return True

            curr = curr.nextHE
            if(curr == first):
                break

        return False


    # Iterate over the verts that make up this face
    def adjVerts(self):

        # Iterate through the adjacent faces
        first = self.anyHalfEdge
        curr = self.anyHalfEdge
        while True:

            # Yield a vertex
            yield curr.vertex

            curr = curr.nextHE
            if(curr == first):
                break

    # Iterate over the halfedges that make up this face
    def adjHalfEdges(self):

        # Iterate through the adjacent faces
        first = self.anyHalfEdge
        curr = self.anyHalfEdge
        while True:

            # Yield a halfedge
            yield curr

            curr = curr.nextHE
            if(curr == first):
                break

    # Iterate over the edges that make up this face
    def adjEdges(self):

        # Iterate through the adjacent faces
        first = self.anyHalfEdge
        curr = self.anyHalfEdge
        while True:

            # Yield an edge
            yield curr.edge

            curr = curr.nextHE
            if(curr == first):
                break


    @property
    def normal(self):
        """The normal vector for this face"""
        if 'normal' in self._cache: return self._cache['normal']

        v = list(self.adjVerts())
        n = normalize(cross(v[1].pos - v[0].pos, v[2].pos - v[0].pos))

        if self.staticGeometry: self._cache['normal'] = n
        return n


    @property
    def area(self):
        """The area of this face"""
        if 'area' in self._cache: return self._cache['area']

        v = list(self.adjVerts())
        a = 0.5 * norm(cross(v[1].pos - v[0].pos, v[2].pos - v[0].pos))

        if self.staticGeometry: self._cache['area'] = a
        return a

class Edge(object):


    def __init__(self, staticGeometry=False):
        ### Members
        anyHalfEdge = None      # Either of the halfedges (if this is a boundary edge,
                                # guaranteed to be the real one)

        self._cache = dict()
        self.staticGeometry = staticGeometry


    # Return a boolean indicating whether this is on the boundary of the mesh
    def isBoundary(self):
        return self.anyHalfEdge.isBoundary()

    # Return true if the edge can be flipped, meaning:
    #   - The edge is not on the boundary
    #   - Neither of the verts that neighbor the edge has degree <= 3
    def canFlip(self):

        if self.isBoundary():
            return False

        # Can only flip if both vertices have degree > 3
        v1 = self.anyHalfEdge.vertex
        v2 = self.anyHalfEdge.twinHE.vertex

        return v1.degree() > 3 and v2.degree() > 3


    # Flip this edge
    # Does nothing if canFlip() returns false
    def flip(self):

        if self.staticGeometry:
            raise ValueError("ERROR: Cannot flip edge with static geometry")

        if not self.canFlip():
            return

        # Note: This does a complete reassignment of references, which will likely include
        # changing things that don't technically need to be changed (like anyHalfEdge refs
        # that weren't actually invalidated). This is done for simplicity and conciseness.


        # Get references to the relevant objects
        h1 = self.anyHalfEdge
        h2 = h1.twinHE
        h1n = h1.nextHE
        h1nn = h1n.nextHE
        h2n = h2.nextHE
        h2nn = h2n.nextHE
        e = h1.edge
        f1 = h1.face
        f2 = h2.face
        Vold1 = h1.vertex
        Vold2 = h2.vertex
        Vnew1 = h1.nextHE.vertex
        Vnew2 = h2.nextHE.vertex


        ## Re-assign pointers
        # This edge
        self.anyHalfEdge = h1
        h1.vertex = Vnew2
        h2.vertex = Vnew1

        # Lower face HE loop
        h1.nextHE = h2nn
        h2nn.nextHE = h1n
        h1n.nextHE = h1

        # Upper face HE loop
        h2.nextHE = h1nn
        h1nn.nextHE = h2n
        h2n.nextHE = h2

        # Faces
        f1.anyHalfEdge = h1
        f2.anyHalfEdge = h2
        h2nn.face = f1
        h1nn.face = f2

        # Verts
        Vold1.anyHalfEdge = h1n
        Vold2.anyHalfEdge = h2n


    @property
    def cotanLaplace(self):
        """
        Return the cotan-laplace weight of the this edge
        """
        if 'cotanLaplace' in self._cache: return self._cache['cotanLaplace']

        # TODO implement me
        val = 0.0

        if self.staticGeometry: self._cache['cotanLaplace'] = val
        return val
