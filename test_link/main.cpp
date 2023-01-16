#include <iostream>
/*----------------------------------------------------------------------------*/
#include <Triton2/TetgenInterface/TetgenFacade.h>
/*----------------------------------------------------------------------------*/
using namespace gmds;
int main() {
        MeshModel mod = DIM2|N|F|F2N;
        IGMesh boundary_mesh(mod);
        boundary_mesh.newNode(0,0);
        boundary_mesh.newNode(0,1);
        boundary_mesh.newNode(1,0);
        boundary_mesh.newNode(1,1);

        IGMesh vol_mesh(mod);
        TetgenFacade f;
        f.generateTetMesh(boundary_mesh, vol_mesh);

        return 0;
}
