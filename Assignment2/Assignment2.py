# Basic application to load a mesh from file and view it in a window

# Python imports
import sys, os
import euclid as eu

## Imports from this project
sys.path.append(os.path.join(os.path.dirname(__file__), '..', 'core')) # hack to allow local imports without creaing a module or modifying the path variable
from InputOutput import *
from MeshDisplay import MeshDisplay
from HalfEdgeMesh import *
from Utilities import *

def main():

    # Get the path for the mesh to load, either from the program argument if
    # one was given, or a dialog otherwise
    if(len(sys.argv) > 1):
        filename = sys.argv[1]
    else:
        print("ERROR: No file name specified. Proper syntax is 'python Assignment2.py path/to/your/mesh.obj'.")
        exit()

    # Read in the mesh
    mesh = HalfEdgeMesh(readMesh(filename))

    # Create a viewer object
    winName = 'DDG Assignment2 -- ' + os.path.basename(filename)
    meshDisplay = MeshDisplay(windowTitle=winName, width=400, height=300)
    meshDisplay.setMesh(mesh)

    ###################### BEGIN YOUR CODE
    # implement the body of each of these functions

    @property
    @cacheGeometry
    def faceArea(self):
        """
        Compute the area of a face. Though not directly requested, this will be
        useful when computing face-area weighted normals below.
        This method gets called on a face, so 'self' is a reference to the
        face at which we will compute the area.
        """
        v = list(self.adjacentVerts())
        a = 0.5 * norm(cross(v[1].position - v[0].position, v[2].position - v[0].position))

        return a

    def faceArea2(self):
        """
        use area vector to compute the polygon area
        """
        sum_areavector = [0.0, 0.0, 0.0]

        verts = list(self.adjacentVerts())
        LEN = len(verts)
        for (i, v) in enumerate(verts):
            sum_areavector += 0.5 * cross(verts[i].position, verts[(i+1)%LEN].position)

        return norm(sum_areavector)

    @property
    @cacheGeometry
    def faceNormal(self):
        """
        Compute normal at a face of the mesh. Unlike at vertices, there is one very
        obvious way to do this, since a face uniquely defines a plane.
        This method gets called on a face, so 'self' is a reference to the
        face at which we will compute the normal.
        """
        v = list(self.adjacentVerts())
        n = normalize(cross(v[1].position - v[0].position, v[2].position - v[0].position))

        return n


    @property
    @cacheGeometry
    def vertexNormal_EquallyWeighted(self):
        """
        Compute a vertex normal using the 'equally weighted' method.
        This method gets called on a vertex, so 'self' is a reference to the
        vertex at which we will compute the normal.
        """
        normalSum = np.array([0.0,0.0,0.0])
        for face in self.adjacentFaces():
            normalSum += face.normal * 1.0
        n = normalize(normalSum)

        return n

    @property
    @cacheGeometry
    def vertexNormal_AreaWeighted(self):
        """
        Compute a vertex normal using the 'face area weights' method.
        This method gets called on a vertex, so 'self' is a reference to the
        vertex at which we will compute the normal.
        """
        normalSum = np.array([0.0,0.0,0.0])
        for face in self.adjacentFaces():
            normalSum += face.normal * face.area
        n = normalize(normalSum)

        return n

    @property
    @cacheGeometry
    def vertexNormal_AngleWeighted(self):
        """
        Compute a vertex normal using the 'tip angle weights' method.
        This method gets called on a vertex, so 'self' is a reference to the
        vertex at which we will compute the normal.
        """
        normalSum = np.array([0.0,0.0,0.0])

        for face in self.adjacentFaces():

            v = list(face.adjacentVerts())

            v0 = v1 = v2 = Vertex()

            if v[0].id == self.id:
               v0 = v[0]; v1 = v[1]; v2 = v[2];
            if v[1].id == self.id:
               v0 = v[1]; v1 = v[2]; v2 = v[0];
            if v[2].id == self.id:
               v0 = v[2]; v1 = v[0]; v2 = v[1];

            a = v1.position - v0.position
            b = v2.position - v0.position

            theta = acos(np.dot((a/norm(a)),(b/norm(b))))

            normalSum += face.normal * theta

        n = normalize(normalSum)
        return n


    #@property
    #@cacheGeometry
    def cotan(self):
        """
        Compute the cotangent of the angle opposite a halfedge. This is not
        directly required, but will be useful when computing the mean curvature
        normals below.
        This method gets called on a halfedge, so 'self' is a reference to the
        halfedge at which we will compute the cotangent.
        """
        if self.next.next.next is not self:
            raise ValueError("ERROR: halfedge.cotan() is only well-defined on a triangle")

        if self.isReal:

            # Relevant vectors
            v0 =  self.next.next.vector
            v1 = -self.next.vector

            # Nifty vector equivalent of cot(theta)
            val = np.dot(v0, v1) / norm(cross(v0, v1))
            return val

        else:
            return 0.0 # placeholder value


    @property
    @cacheGeometry
    def vertexNormal_MeanCurvature(self):
        """
        Compute a vertex normal using the 'mean curvature' method.
        Be sure to understand the relationship between this method and the
        area gradient method.
        This method gets called on a vertex, so 'self' is a reference to the
        vertex at which we will compute the normal.
        """
        # the vertex normals we get from the mean curvature vector are precisely 
        # the same as the ones we get from the area gradient
        # areaGrad(p_i) = 0.5 * SUM((cot a_j + cot b_j)(p_j - p_i))
        sum_normal = [0.0, 0.0, 0.0]

        halfedges = list(self.adjacentHalfEdges_CounterClockwise())
        for he in halfedges:
            sum_normal += (cotan(he.twin) + cotan(he)) * he.vector # (p_j - p_i) = he.vector

        n = normalize(0.5 * sum_normal)
        return -n # But it seems that I should use -n to return the correct normal value

    @property
    @cacheGeometry
    def vertexNormal_SphereInscribed(self):
        """
        Compute a vertex normal using the 'inscribed sphere' method.
        This method gets called on a vertex, so 'self' is a reference to the
        vertex at which we will compute the normal.
        """
        # Ns = 1/c * SUM(e(j) x e(j+1) / (|e(j)|^2 * |e(j+1)|^2))
        sum_normal = [0.0, 0.0, 0.0]

        halfedges = list(self.adjacentHalfEdges_CounterClockwise())
        LEN = len(halfedges)

        for j in range(0, LEN-1): # [0, LEN-1)
            normj  = norm(halfedges[j].vector)
            normj1 = norm(halfedges[j+1].vector)

            sum_normal += cross(halfedges[j].vector, halfedges[j+1].vector) / (normj*normj * normj1*normj1)

        n = normalize(sum_normal)
        return n # But it seems that I should use -n to return the correct normal value



    @property
    @cacheGeometry
    def angleDefect(self):
        """
        Compute the angle defect of a vertex, d(v) (see Assignment 1 Exercise 8).
        This method gets called on a vertex, so 'self' is a reference to the
        vertex at which we will compute the angle defect.
        """
        sum_theta = 0.0
        # acos(np.dot((a/norm(a)),(b/norm(b))))
        for face in self.adjacentFaces():
            
            v = list(face.adjacentVerts())

            v0 = v1 = v2 = Vertex()

            if v[0].id == self.id:
               v0 = v[0]; v1 = v[1]; v2 = v[2];
            if v[1].id == self.id:
               v0 = v[1]; v1 = v[2]; v2 = v[0];
            if v[2].id == self.id:
               v0 = v[2]; v1 = v[0]; v2 = v[1];

            a = v1.position - v0.position
            b = v2.position - v0.position

            theta = acos(np.dot((a/norm(a)),(b/norm(b))))

            sum_theta += theta
	
        return 2.0 * pi - sum_theta


    def totalGaussianCurvature():
        """
        Compute the total Gaussian curvature in the mesh, meaning the sum of Gaussian
        curvature at each vertex.
        Note that you can access the mesh with the 'mesh' variable.
        """
        sum_ = 0.0
        for v in mesh.verts:
            sum_ += v.angleDefect
        return sum_


    def gaussianCurvatureFromGaussBonnet():
        """
        Compute the total Gaussian curvature that the mesh should have, given that the
        Gauss-Bonnet theorem holds (see Assignment 1 Exercise 9).
        Note that you can access the mesh with the 'mesh' variable. The
        mesh includes members like 'mesh.verts' and 'mesh.faces', which are
        sets of the vertices (resp. faces) in the mesh.
        """
        X = len(mesh.verts) - len(mesh.edges) + len(mesh.faces)
        return 2.0 * pi * X


    ###################### END YOUR CODE


    # Set these newly-defined methods as the methods to use in the classes
    Face.normal = faceNormal
    Face.area = faceArea
    Vertex.normal = vertexNormal_AreaWeighted
    Vertex.vertexNormal_EquallyWeighted = vertexNormal_EquallyWeighted
    Vertex.vertexNormal_AreaWeighted = vertexNormal_AreaWeighted
    Vertex.vertexNormal_AngleWeighted = vertexNormal_AngleWeighted
    Vertex.vertexNormal_MeanCurvature = vertexNormal_MeanCurvature
    Vertex.vertexNormal_SphereInscribed = vertexNormal_SphereInscribed
    Vertex.angleDefect = angleDefect
    HalfEdge.cotan = cotan


    ## Functions which will be called by keypresses to visualize these definitions

    def toggleFaceVectors():
        print("\nToggling vertex vector display")
        if toggleFaceVectors.val:
            toggleFaceVectors.val = False
            meshDisplay.setVectors(None)
        else:
            toggleFaceVectors.val = True
            meshDisplay.setVectors('normal', vectorDefinedAt='face')
        meshDisplay.generateVectorData()
    toggleFaceVectors.val = False # ridiculous Python scoping hack
    meshDisplay.registerKeyCallback('1', toggleFaceVectors, docstring="Toggle drawing face normal vectors")


    def toggleVertexVectors():
        print("\nToggling vertex vector display")
        if toggleVertexVectors.val:
            toggleVertexVectors.val = False
            meshDisplay.setVectors(None)
        else:
            toggleVertexVectors.val = True
            meshDisplay.setVectors('normal', vectorDefinedAt='vertex')
        meshDisplay.generateVectorData()
    toggleVertexVectors.val = False # ridiculous Python scoping hack
    meshDisplay.registerKeyCallback('2', toggleVertexVectors, docstring="Toggle drawing vertex normal vectors")


    def toggleDefect():
        print("\nToggling angle defect display")
        if toggleDefect.val:
            toggleDefect.val = False
            meshDisplay.setShapeColorToDefault()
        else:
            toggleDefect.val = True
            meshDisplay.setShapeColorFromScalar("angleDefect", cmapName="seismic",vMinMax=[-pi/8,pi/8])
        meshDisplay.generateFaceData()
    toggleDefect.val = False
    meshDisplay.registerKeyCallback('3', toggleDefect, docstring="Toggle drawing angle defect coloring")


    def useEquallyWeightedNormals():
        mesh.staticGeometry = False
        print("\nUsing equally-weighted normals")
        Vertex.normal = vertexNormal_EquallyWeighted
        mesh.staticGeometry = True
        meshDisplay.generateAllMeshValues()
    meshDisplay.registerKeyCallback('4', useEquallyWeightedNormals, docstring="Use equally-weighted normal computation")

    def useAreaWeightedNormals():
        mesh.staticGeometry = False
        print("\nUsing area-weighted normals")
        Vertex.normal = vertexNormal_AreaWeighted
        mesh.staticGeometry = True
        meshDisplay.generateAllMeshValues()
    meshDisplay.registerKeyCallback('5', useAreaWeightedNormals, docstring="Use area-weighted normal computation")

    def useAngleWeightedNormals():
        mesh.staticGeometry = False
        print("\nUsing angle-weighted normals")
        Vertex.normal = vertexNormal_AngleWeighted
        mesh.staticGeometry = True
        meshDisplay.generateAllMeshValues()
    meshDisplay.registerKeyCallback('6', useAngleWeightedNormals, docstring="Use angle-weighted normal computation")

    def useMeanCurvatureNormals():
        mesh.staticGeometry = False
        print("\nUsing mean curvature normals")
        Vertex.normal = vertexNormal_MeanCurvature
        mesh.staticGeometry = True
        meshDisplay.generateAllMeshValues()
    meshDisplay.registerKeyCallback('7', useMeanCurvatureNormals, docstring="Use mean curvature normal computation")

    def useSphereInscribedNormals():
        mesh.staticGeometry = False
        print("\nUsing sphere-inscribed normals")
        Vertex.normal = vertexNormal_SphereInscribed
        mesh.staticGeometry = True
        meshDisplay.generateAllMeshValues()
    meshDisplay.registerKeyCallback('8', useSphereInscribedNormals, docstring="Use sphere-inscribed normal computation")

    def computeDiscreteGaussBonnet():
        print("\nComputing total curvature:")
        computed = totalGaussianCurvature()
        predicted = gaussianCurvatureFromGaussBonnet()
        print("   Total computed curvature: " + str(computed))
        print("   Predicted value from Gauss-Bonnet is: " + str(predicted))
        print("   Error is: " + str(abs(computed - predicted)))
    meshDisplay.registerKeyCallback('z', computeDiscreteGaussBonnet, docstring="Compute total curvature")

    def deformShape():
        print("\nDeforming shape")
        mesh.staticGeometry = False

        # Get the center and scale of the shape
        center = meshDisplay.dataCenter
        scale = meshDisplay.scaleFactor

        # Rotate according to swirly function
        ax = eu.Vector3(-1.0,.75,0.5)
        for v in mesh.verts:
            vec = v.position - center
            theta = 0.8 * norm(vec) / scale
            newVec = np.array(eu.Vector3(*vec).rotate_around(ax, theta))
            v.position = center + newVec


        mesh.staticGeometry = True
        meshDisplay.generateAllMeshValues()

    meshDisplay.registerKeyCallback('x', deformShape, docstring="Apply a swirly deformation to the shape")



    ## Register pick functions that output useful information on click
    def pickVert(vert):
        print("   Position:" + printVec3(vert.position))
        print("   Angle defect: {:.5f}".format(vert.angleDefect))
        print("   Normal (equally weighted): " + printVec3(vert.vertexNormal_EquallyWeighted))
        print("   Normal (area weighted):    " + printVec3(vert.vertexNormal_AreaWeighted))
        print("   Normal (angle weighted):   " + printVec3(vert.vertexNormal_AngleWeighted))
        print("   Normal (sphere-inscribed): " + printVec3(vert.vertexNormal_SphereInscribed))
        print("   Normal (mean curvature):   " + printVec3(vert.vertexNormal_MeanCurvature))
    meshDisplay.pickVertexCallback = pickVert

    def pickFace(face):
        print("   Face area : {:.5f}".format(face.area))
        print("   Face area2: {:.5f}".format(faceArea2(face)))
        print("   Normal: " + printVec3(face.normal))
        print("   Vertex positions: ")
        for (i, vert) in enumerate(face.adjacentVerts()):
            print("     v{}: {}".format((i+1),printVec3(vert.position)))
    meshDisplay.pickFaceCallback = pickFace


    # Start the viewer running
    meshDisplay.startMainLoop()


if __name__ == "__main__": main()
