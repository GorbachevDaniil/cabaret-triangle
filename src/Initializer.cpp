#include "Initializer.hpp"

void Initializer::initialize(Mesh &mesh) {
    for (unsigned long i = 0; i < mesh.nodes.size(); i++) {
        Data *data = &mesh.nodes[i].data;
        data->u = Vector(1, 1);
        if (data->coords.x < -0.5 && data->coords.x > -0.9 &&
                data->coords.y < -0.5 && data->coords.y > -0.9) {
            data->phi0 = 1;
        } else {
            data->phi0 = 0;
        }
    }
}
