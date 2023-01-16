#include <iostream>
/*----------------------------------------------------------------------------*/
#include <Triton2/NetgenInterface/NetgenFacade.h>
/*----------------------------------------------------------------------------*/
using namespace gmds;
int main() {
        MeshModel mod = DIM2|N|F|F2N;
        IGMesh mesh(mod);
        mesh.newNode(0,0);
        mesh.newNode(0,1);
        mesh.newNode(1,0);
        mesh.newNode(1,1);

        NetgenFacade f;
        f.createGMDSMesh(mesh);

        return 0;
}
