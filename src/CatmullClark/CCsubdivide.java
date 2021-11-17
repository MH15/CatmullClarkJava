package CatmullClark;

import osu.halfEdgeMesh.*;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.List;
import java.util.stream.Collectors;

public class CCsubdivide {

    public static void main(String[] args) {
        String inputFilename = "data/cubeQ.off";
//        String inputFilename = "data/quad5A.off";
//        String inputFilename = "data/torus-isov0.1-quads.off";
        var oldMesh = new HalfEdgeMeshA();
        var reader = new OffFileReaderA();
        reader.OpenAndReadFile(inputFilename, oldMesh);

        int steps = 3;

        if (steps > 0) {

            var newMesh = subdivide(oldMesh);
            for (int i = 1; i < steps; i++) {
                newMesh = subdivide(newMesh);

            }
            var writer = new OffFileWriterA();
            writer.OpenAndWriteFile("data/test.off", newMesh);
        } else {
            var writer = new OffFileWriterA();
            writer.OpenAndWriteFile("data/test.off", oldMesh);
        }


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

                if (createdOriginals.containsKey(oldVertexIndex)) {
                    cellVertIndices.get(cellIndex).add(createdOriginals.get(oldVertexIndex));
                } else {
                    if (!edge.IsBoundary()) {
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

                        float[] newPoint = sumCoords(weightedOriginalCoords, weightedFacePoints, weightedEdgePoints);

                        int newFromVertIndex = createVertexAtPoint(newMesh, newPoint);

                        cellVertIndices.get(cellIndex).add(newFromVertIndex);
                        createdOriginals.put(oldVertexIndex, newFromVertIndex);
                    } else {
                        int n = fromVertex.NumHalfEdgesFrom();

                        float[][] eArray = new float[n][3];
                        int eArrayCount = n;

                        float[][] fArray = new float[n][3];
                        for (int i = 0; i < n; i++) {
                            var e = fromVertex.KthHalfEdgeFrom(i);
                            if (e.IsBoundary()) {
                                eArrayCount--;
                            }

                            fArray[i] = facePointEdges.get(e.Index()).coord;
                            eArray[i] = averageCoords(fromVertex.coord, e.ToVertex().coord);
                        }

                        float[][] eArrayWeCareAbout = Arrays.copyOfRange(eArray, 0, eArrayCount);

                        var F = divCoord(sumCoords(fArray), eArrayCount);
                        var E = divCoord(sumCoords(eArrayWeCareAbout), eArrayCount);
                        System.out.println(eArrayCount);


                        float m1 = (n - 3) / (float) (n);
                        float m2 = 1 / (float) (eArrayCount);
                        float m3 = 2 / (float) (eArrayCount);

                        float[] weightedOriginalCoords = mulCoord(fromVertex.coord, m1);
                        float[] weightedFacePoints = mulCoord(F, m2);
                        float[] weightedEdgePoints = mulCoord(E, m3);

                        float[] newPoint = averageCoords(fromVertex.coord, E);

                        int newFromVertIndex = createVertexAtPoint(newMesh, newPoint);

                        System.out.println(Arrays.toString(E));
                        cellVertIndices.get(cellIndex).add(newFromVertIndex);
                        createdOriginals.put(oldVertexIndex, newFromVertIndex);
                    }
                }


                int oppositeIndex = edge.NextHalfEdgeAroundEdge().Index();
                if (!splitEdges.containsKey(edgeIndex)) {
                    float startCoords[] = fromVertex.coord;
                    float endCoords[] = edge.ToVertex().coord;

                    var facePointA = facePointEdges.get(edgeIndex);
                    var facePointB = facePointEdges.get(oppositeIndex);

                    float newCoords[] = averageCoords(startCoords, endCoords);

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

        for (
                var cellIndex : oldMesh.CellIndices()) {
            var cell = oldMesh.Cell(cellIndex);

            var verts = cellVertIndices.get(cellIndex);

            int numVerts = verts.size();

            int offset = numVerts * 2 - 2;
            for (int i = 1; i < verts.size(); i += 2) {
//                int v0 = verts.get((i + 14) % numVerts);
                int v0 = verts.get((i + offset) % numVerts);
                int v1 = verts.get((i - 1) % numVerts);
                int v2 = verts.get((i + 0) % numVerts);
                int v3 = cellIndexToFacePointIndex.get(cellIndex);
                createCell(newMesh, v0, v1, v2, v3);

            }
        }
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
