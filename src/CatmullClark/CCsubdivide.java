package CatmullClark;

import osu.halfEdgeMesh.*;
import picocli.CommandLine;

import java.util.*;
import java.util.concurrent.Callable;
import java.util.stream.Collectors;

@CommandLine.Command(name = "CCsubdivide", mixinStandardHelpOptions = true, version = "0.1",
        description = "Subdivide an OFF file using the Catmull-Clark algorithm.")
public class CCsubdivide implements Callable<Integer> {

    @CommandLine.Parameters(index = "0", description = "Input file in .off format.", paramLabel = "<input .off file>")
    private String filename;


    @CommandLine.Option(names = {"-n", "-num_iter"}, description = "integer 1 or greater")
    private int n;

    @Override
    public Integer call() {
        String outputFilename = "out.off";

        var oldMesh = new HalfEdgeMeshA();
        var reader = new OffFileReaderA();
        reader.OpenAndReadFile(filename, oldMesh);


        if (n > 0) {
            System.out.println("Iteration 0");
            var newMesh = subdivide(oldMesh);
            checkMesh(newMesh);

            for (int i = 1; i < n; i++) {
                System.out.println("Iteration " + i);
                newMesh = subdivide(newMesh);
                checkMesh(newMesh);
            }

            var writer = new OffFileWriterA();
            writer.OpenAndWriteFile(outputFilename, newMesh);
        } else {
            var writer = new OffFileWriterA();
            writer.OpenAndWriteFile(outputFilename, oldMesh);
        }

        return 0;
    }

    public static void main(String[] args) {
        // Parse command line options with Picocli
        int exitCode = new CommandLine(new CCsubdivide()).execute(args);
        System.exit(exitCode);
    }

    static HalfEdgeMeshA subdivide(HalfEdgeMeshA oldMesh) {
        var newMesh = new HalfEdgeMeshA();

        // Edge index to Vertex
        var splitEdges = new HashMap<Integer, VertexBase>();

        // Edge index to Vertex
        var facePointEdges = new HashMap<Integer, VertexBase>();

        // Cell index to face point
        var cellIndexToFacePointIndex = new HashMap<Integer, Integer>();

        // Create the face points for each cell
        for (var cellIndex : oldMesh.CellIndices()) {
            var cell = oldMesh.Cell(cellIndex);

            var edges = getEdgesOfCell(cell);
            int facePointIndex = addFacePoint(newMesh, edges);
            for (var edge : edges) {
                facePointEdges.put(edge.Index(), newMesh.Vertex(facePointIndex));
                cellIndexToFacePointIndex.put(cellIndex, facePointIndex);
            }

        }


        // Old vertex index, new vertex index
        var createdOriginals = new HashMap<Integer, Integer>();

        // Cell index to Vertex Indices
        var cellVertIndices = new HashMap<Integer, List<Integer>>();


        for (var cellIndex : oldMesh.CellIndices()) {
            var cell = oldMesh.Cell(cellIndex);
            cellVertIndices.put(cellIndex, new ArrayList<>());

            var edges = getEdgesOfCell(cell);
            for (var edge : edges) {


                int edgeIndex = edge.Index();
                VertexBase fromVertex = edge.FromVertex();

                int oldVertexIndex = fromVertex.Index();

                // If vertex is on the boundary, don't move it
//                if (edge.IsBoundary() || edge.PrevHalfEdgeInCell().IsBoundary()) {
                if (fromVertex.KthHalfEdgeFrom(0).IsBoundary() && !createdOriginals.containsKey(oldVertexIndex)) {
                    int newFromVertIndex = createVertexAtPoint(newMesh, fromVertex.coord);

                    cellVertIndices.get(cellIndex).add(newFromVertIndex);
                    createdOriginals.put(oldVertexIndex, newFromVertIndex);

                }
                // If the point has already been created, retrieve it
                else if (createdOriginals.containsKey(oldVertexIndex)) {
                    cellVertIndices.get(cellIndex).add(createdOriginals.get(oldVertexIndex));
                }
                // Otherwise, need to calculate new point position
                else {
                    float[] newPoint = findNewVertexPoints(fromVertex, facePointEdges);
                    int newFromVertIndex = createVertexAtPoint(newMesh, newPoint);

                    cellVertIndices.get(cellIndex).add(newFromVertIndex);
                    createdOriginals.put(oldVertexIndex, newFromVertIndex);
                }


                int oppositeIndex = edge.NextHalfEdgeAroundEdge().Index();
                if (!splitEdges.containsKey(edgeIndex)) {
                    float[] startCoords = fromVertex.coord;
                    float[] endCoords = edge.ToVertex().coord;

                    var facePointA = facePointEdges.get(edgeIndex);
                    var facePointB = facePointEdges.get(oppositeIndex);

                    float[] newCoords = averageCoords(startCoords, endCoords);

                    // "for the edges that are on the border of a hole, the edge point is just the middle of the edge."
                    if (!edge.IsBoundary()) {

                        if (facePointB != null) {
                            newCoords = averageCoords(startCoords, endCoords, facePointA.coord, facePointB.coord);
                        } else {
                            newCoords = averageCoords(startCoords, endCoords, facePointA.coord);
                        }
                    }

                    int newVertIndex = createVertexAtPoint(newMesh, newCoords);
                    cellVertIndices.get(cellIndex).add(newVertIndex);
                    var vertex = newMesh.Vertex(newVertIndex);
                    splitEdges.put(edgeIndex, vertex);
                    splitEdges.put(edge.NextHalfEdgeAroundEdge().Index(), vertex);
                } else {
                    var vertex = splitEdges.get(edgeIndex);
                    cellVertIndices.get(cellIndex).add(vertex.Index());
                }
            }


        }

        // Add the faces
        for (var cellIndex : oldMesh.CellIndices()) {
            var verts = cellVertIndices.get(cellIndex);
            int offset = verts.size() * 2 - 2;
            for (int i = 1; i < verts.size(); i += 2) {
                int v0 = verts.get((i + offset) % verts.size());
                int v1 = verts.get((i - 1) % verts.size());
                int v2 = verts.get((i) % verts.size());
                int v3 = cellIndexToFacePointIndex.get(cellIndex);
                createCell(newMesh, v0, v1, v2, v3);
            }
        }
        return newMesh;
    }

    /**
     * Validate created mesh
     *
     * @param mesh mesh to check
     */
    static void checkMesh(HalfEdgeMeshA mesh) {
        var manifold_info = mesh.CheckManifold();

        var error_info = mesh.CheckAll();
        if (error_info.FlagError()) {
            System.err.println("Error detected in mesh data structure.");
            if (error_info.Message() != null && !(error_info.Message().equals(""))) {
                System.err.println(error_info.Message());
            }
            System.exit(-1);
        }

        if (!manifold_info.FlagManifold()) {
            if (!manifold_info.FlagManifoldVertices()) {
                int iv = manifold_info.VertexIndex();
                System.err.println
                        ("Warning: Non manifold vertex " +
                                iv + ".");
            }

            if (!manifold_info.FlagManifoldEdges()) {
                int ihalf_edge = manifold_info.HalfEdgeIndex();
                osu.halfEdgeMesh.HalfEdgeBase half_edge = mesh.HalfEdge(ihalf_edge);
                System.err.println
                        ("Warning: Non manifold edge (" +
                                half_edge.EndpointsStr(",") + ").");
            }

        }
    }


    /**
     * Find the vertex locations of existing points, per the Catmull-Clark algorithm.
     *
     * @param fromVertex     current point on the old mesh
     * @param facePointEdges mapping edge index to vertex index
     * @return new coordinate position
     */
    static float[] findNewVertexPoints(VertexBase fromVertex, Map<Integer, VertexBase> facePointEdges) {
        int n = fromVertex.NumHalfEdgesFrom();

        float[][] rArray = new float[n][3];

        float[][] fArray = new float[n][3];
        for (int i = 0; i < n; i++) {
            var e = fromVertex.KthHalfEdgeFrom(i);
            fArray[i] = facePointEdges.get(e.Index()).coord;
            rArray[i] = averageCoords(fromVertex.coord, e.ToVertex().coord);
        }

        var F = averageCoords(fArray);
        var R = averageCoords(rArray);


        float m1 = (n - 3) / (float) (n);
        float m2 = 1 / (float) (n);
        float m3 = 2 / (float) (n);

        float[] weightedOriginalCoords = mulCoord(fromVertex.coord, m1);
        float[] weightedFacePoints = mulCoord(F, m2);
        float[] weightedEdgePoints = mulCoord(R, m3);

        return sumCoords(weightedOriginalCoords, weightedFacePoints, weightedEdgePoints);
    }

    /**
     * Add a new cell to the mesh
     *
     * @param mesh   mesh the cell is added to
     * @param vIndex one or more vertices to form the cell from
     */
    static void createCell(HalfEdgeMeshA mesh, int... vIndex) {
        int cellIndex = mesh.MaxCellIndex() + 1;

        ArrayList<Integer> list = (ArrayList<Integer>) Arrays.stream(vIndex).boxed().collect(Collectors.toList());
        try {
            mesh.AddCell(cellIndex, list);
        } catch (Exception e) {
            e.printStackTrace();
        }
    }


    /**
     * Return the coordinates multiplied by a scalar
     *
     * @param a      coordinate
     * @param scalar float
     * @return new coordinates
     */
    static float[] mulCoord(float[] a, float scalar) {
        float[] result = {0, 0, 0};
        for (int j = 0; j < 3; j++) {
            result[j] = a[j] * scalar;
        }
        return result;
    }


    /**
     * Sum a list of float[3]
     *
     * @param a one or more coordinates
     * @return sum value
     */
    static float[] sumCoords(float[]... a) {
        float[] result = {0, 0, 0};
        for (float[] floats : a) {
            for (int j = 0; j < 3; j++) {
                result[j] += floats[j];
            }
        }
        return result;
    }

    /**
     * Average a list of float[3]
     *
     * @param a one or more coordinates
     * @return averaged value
     */
    static float[] averageCoords(float[]... a) {
        float[] result = sumCoords(a);
        result[0] /= a.length;
        result[1] /= a.length;
        result[2] /= a.length;
        return result;
    }

    /**
     * Walk around the half-edges in a cell and build
     * an ordered list of all edges in said cell.
     *
     * @param cell the face
     * @return list of edges surrounding the face
     */
    static List<HalfEdgeBase> getEdgesOfCell(CellBase cell) {
        // Gather all edges of the cell
        ArrayList<HalfEdgeBase> edgesOfCell = new ArrayList<>();
        var halfEdge = cell.HalfEdge();
        edgesOfCell.add(halfEdge);
        var nextEdge = halfEdge.NextHalfEdgeInCell();
        while (!nextEdge.equals(halfEdge)) {
            edgesOfCell.add(nextEdge);
            nextEdge = nextEdge.NextHalfEdgeInCell();
        }

        return edgesOfCell;
    }

    /**
     * Add a vertex at the average of a set of points.
     *
     * @param mesh  mesh to add vertex to
     * @param edges edges to average
     * @return index of new point
     */
    static int addFacePoint(HalfEdgeMeshA mesh, List<HalfEdgeBase> edges) {
        int numEdges = edges.size();

        float[] coords = {0, 0, 0};

        for (var edge : edges) {
            var vert = edge.FromVertex();
            coords[0] += vert.coord[0];
            coords[1] += vert.coord[1];
            coords[2] += vert.coord[2];
        }
        coords[0] /= numEdges;
        coords[1] /= numEdges;
        coords[2] /= numEdges;
        return createVertexAtPoint(mesh, coords);
    }

    /**
     * Create a new vertex in the mesh.
     *
     * @param mesh   mesh to add vertex to
     * @param coords where to add the vertex
     * @return index of newly created vertex
     */
    static int createVertexAtPoint(HalfEdgeMeshA mesh, float[] coords) {
        try {
            mesh.AddVertices(0);
            int vertIndex = mesh.MaxVertexIndex() + 1;
            mesh.SetCoord(vertIndex, coords);
            return vertIndex;
        } catch (Exception e) {
            e.printStackTrace();
        }
        return -1;
    }


}
