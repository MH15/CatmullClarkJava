package CatmullClark;

import osu.halfEdgeMesh.*;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.List;
import java.util.stream.Collectors;

public class CCsubdivide {

    public static void main(String[] args) {
        String inputFilename = "data/cubeT.off";
//        String inputFilename = "data/torus-isov0.1-quads.off";
        var oldMesh = new HalfEdgeMeshA();
        var reader = new OffFileReaderA();
        reader.OpenAndReadFile(inputFilename, oldMesh);


        var newMesh = subdivide(oldMesh);
//        newMesh = subdivide(newMesh);
//        newMesh = subdivide(newMesh);

        var writer = new OffFileWriterA();
        writer.OpenAndWriteFile("data/test.off", newMesh);


    }

    static HalfEdgeMeshA subdivide(HalfEdgeMeshA oldMesh) {
        var newMesh = new HalfEdgeMeshA();

        // Edge index to Vertex
        var splitEdges = new HashMap<Integer, VertexBase>();

        // Edge index to Vertex
        var facePointEdges = new HashMap<Integer, VertexBase>();

        // Cell index to face point
        var cellIndexToFacePointIndex = new HashMap<Integer, Integer>();

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
//                System.out.println("searching for " + vertexIndex);

//                int expectedVertexIndex =

                if (createdOriginals.containsKey(oldVertexIndex)) {
                    cellVertIndices.get(cellIndex).add(createdOriginals.get(oldVertexIndex));
                } else {
                    int n = fromVertex.NumHalfEdgesFrom();
                    var e0 = fromVertex.KthHalfEdgeFrom(0);
                    var e1 = fromVertex.KthHalfEdgeFrom(1);
                    var e2 = fromVertex.KthHalfEdgeFrom(2);

                    var f0 = facePointEdges.get(e0.Index());
                    var f1 = facePointEdges.get(e1.Index());
                    var f2 = facePointEdges.get(e2.Index());

                    var F = divCoord(sumCoords(f0.coord, f1.coord, f2.coord), 3);

                    var s0 = averageCoords(fromVertex.coord, e0.ToVertex().coord);
                    var s1 = averageCoords(fromVertex.coord, e1.ToVertex().coord);
                    var s2 = averageCoords(fromVertex.coord, e2.ToVertex().coord);

                    var R = sumCoords(s0, s1, s2);
                    var R2 = divCoord(mulCoord(R, 2), 3);

                    var FplusR2 = sumCoords(F, R2);

                    float p[] = divCoord(FplusR2, 3);


//                    int newFromVertIndex = createVertexAtPoint(newMesh, p);
                    int newFromVertIndex = createVertexAtPoint(newMesh, fromVertex.coord);
                    cellVertIndices.get(cellIndex).add(newFromVertIndex);
                    createdOriginals.put(oldVertexIndex, newFromVertIndex);
                }


                int oppositeIndex = edge.NextHalfEdgeAroundEdge().Index();
                if (!splitEdges.containsKey(edgeIndex)) {
                    float startCoords[] = fromVertex.coord;
                    float endCoords[] = edge.ToVertex().coord;

                    var facePointA = facePointEdges.get(edgeIndex);
                    var facePointB = facePointEdges.get(oppositeIndex);

                    float newCoords[] = averageCoords(startCoords, endCoords);

                    if (facePointB != null) {
                        newCoords = averageCoords(startCoords, endCoords, facePointA.coord, facePointB.coord);
                    } else {
                        newCoords = averageCoords(startCoords, endCoords, facePointA.coord);
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

//            System.out.println(cellVertIndices);

            int numVerts = cellVertIndices.size();

        }

        for (var cellIndex : oldMesh.CellIndices()) {
            var cell = oldMesh.Cell(cellIndex);

            var verts = cellVertIndices.get(cellIndex);

            int numVerts = verts.size();
            System.out.println("cellVertIndices size: " + numVerts);

            int offset = numVerts * 2 - 2;
            for (int i = 1; i < verts.size(); i += 2) {
//                int v0 = verts.get((i + 14) % numVerts);
                int v0 = verts.get((i + offset) % numVerts);
                int v1 = verts.get((i - 1) % numVerts);
                int v2 = verts.get((i + 0) % numVerts);
                int v3 = cellIndexToFacePointIndex.get(cellIndex);
                System.out.printf("%d -> [%d,%d,%d,%d]%n", i, v0, v1, v2, v3);
                createCell(newMesh, v0, v1, v2, v3);

            }
//            break;
        }
        System.out.println("create");
//        createCell(newMesh, 13, 6, 7, 0);
//        createCell(newMesh, 7, 8, 9, 0);
        return newMesh;
    }

    static void createCell(HalfEdgeMeshA mesh, int... vIndex) {
        int cellIndex = mesh.MaxCellIndex() + 1;

        ArrayList<Integer> list = (ArrayList<Integer>) Arrays.stream(vIndex).boxed().collect(Collectors.toList());
        try {
            mesh.AddCell(cellIndex, list);
        } catch (Exception e) {
            e.printStackTrace();
        }
    }


    static float[] mulCoord(float[] a, float scalar) {
        float[] result = {0, 0, 0};
        for (int j = 0; j < 3; j++) {
            result[j] = a[j] * scalar;
        }
        return result;
    }

    static float[] divCoord(float[] a, float scalar) {
        float[] result = {0, 0, 0};
        for (int j = 0; j < 3; j++) {
            result[j] = a[j] / scalar;
        }
        return result;
    }

    static float[] sumCoords(float[]... a) {
        float[] result = {0, 0, 0};
        for (float[] floats : a) {
            for (int j = 0; j < 3; j++) {
                result[j] += floats[j];
            }
        }
        return result;
    }

    static float[] averageCoords(float[]... a) {
        float[] result = sumCoords(a);
        result[0] /= a.length;
        result[1] /= a.length;
        result[2] /= a.length;
        return result;
    }

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
