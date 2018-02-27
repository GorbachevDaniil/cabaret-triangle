#include "Initializer.hpp"

void Initializer::initialize(Mesh &mesh) {
    for (unsigned long i = 0; i < mesh.nodes.size(); i++) {
        Data *data = &mesh.nodes[i].data;
        data->u.set(1, 1);
        if (data->coords.x < 0.2 && data->coords.x > -0.2 &&
                data->coords.y < 0.2 && data->coords.y > -0.2) {
            data->phi0 = 1;
        } else {
            data->phi0 = 0;
        }
    }
}
