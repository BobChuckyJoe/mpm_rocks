#include <iostream>

#include <Eigen/Dense>

#include "config.hpp"
#include "equations.hpp"

using Eigen::Vector2i;

using namespace std;

struct Particle {
    Vector2d position, velocity;
    double mass, density;
    Matrix2d fe, fp;
};

struct GridCell {
    double mass;
    Vector2d prev_velocity, velocity, force;
};

/*
Given a particle, give the corresponding grid cell it is contained in
*/
Vector2i to_grid(Vector2d position) {
    auto to_ret = (position.array() / GRID_SPACING).floor();
    return Vector2i{(int) to_ret[0], (int) to_ret[1]};
}

double particle_to_grid (Vector2i grid_ind, Vector2d particle_position) {
    return weighting_func(1.0/ GRID_SPACING * (particle_position[0] - grid_ind[0])) *
           weighting_func(1.0/ GRID_SPACING * (particle_position[1] - grid_ind[1]));
}

int main() {
    auto GRID = new GridCell[GRID_LENGTH][GRID_LENGTH];
    auto PARTICLES = new Particle[N_PARTICLES];

    PARTICLES[0].position = Vector2d{3.0, 3.0};
    PARTICLES[1].position = Vector2d{3.0, 4.5};
    PARTICLES[2].position = Vector2d{4.5, 3.0};
    PARTICLES[3].position = Vector2d{4.5, 4.5};
    for (int ind = 0; ind < N_PARTICLES; ind++) {
        PARTICLES[ind].mass = 1.0;
        PARTICLES[ind].density = 1.0;
        PARTICLES[ind].velocity = Vector2d{0.0, 0.0};
        PARTICLES[ind].fe = Matrix2d::Identity();
        PARTICLES[ind].fp = Matrix2d::Identity();
    }

    for (int iteration_num = 0; iteration_num < N_ITERATIONS; iteration_num++) {
        for (int i = 0; i < GRID_LENGTH; i++) {
            for (int j = 0; j < GRID_LENGTH; j++) {
                GRID[i][j].mass = 0.0;
                GRID[i][j].prev_velocity = GRID[i][j].velocity;
                GRID[i][j].velocity = Vector2d{0.0, 0.0};
                GRID[i][j].force = Vector2d{0.0, 0.0};
            }
        }

        // 1
        for (int particle_ind = 0; particle_ind < N_PARTICLES; particle_ind++) {
            auto inds = to_grid(PARTICLES[particle_ind].position);
            for (int dx = -2; dx < 3; dx++) {
                for (int dy = -2; dy < 3; dy++) {
                    auto x = inds[0] + dx;
                    auto y = inds[1] + dy;
                    if (x < 0 || x >= GRID_LENGTH || y < 0 || y >= GRID_LENGTH) {
                        continue;
                    }
                    auto weight = particle_to_grid(Vector2i{x, y}, PARTICLES[particle_ind].position);
                    GRID[x][y].mass += weight * PARTICLES[particle_ind].mass;
                }
            }
        }
        for (int particle_ind = 0; particle_ind < N_PARTICLES; particle_ind++) {\
            auto p = PARTICLES[particle_ind];
            auto inds = to_grid(PARTICLES[particle_ind].position);
            for (int dx = -2; dx < 3; dx++) {
                for (int dy = -2; dy < 3; dy++) {
                    auto x = inds[0] + dx;
                    auto y = inds[1] + dy;
                    if (x < 0 || x >= GRID_LENGTH || y < 0 || y >= GRID_LENGTH) {
                        continue;
                    }
                    auto weight = particle_to_grid(Vector2i{x, y}, PARTICLES[particle_ind].position);
                    auto velocity = PARTICLES[particle_ind].velocity * weight * p.mass / GRID[x][y].mass;
                    GRID[x][y].velocity += velocity;
                }
            }
        }
        //2
        if (iteration_num == 0) {
            for (int particle_ind = 0; particle_ind < N_PARTICLES; particle_ind++) {
                for (int i = 0; i < GRID_LENGTH; i++) {
                    for (int j = 0; j < GRID_LENGTH; j++) {
                        PARTICLES[particle_ind].density += GRID[i][j].mass * particle_to_grid(Vector2i{i, j}, PARTICLES[particle_ind].position);
                    }
                }    
            }
        }
        //3
        for (int i =0; i < N_PARTICLES; i++) {
            auto p = PARTICLES[i];
            auto inds = to_grid(p.position);
            for (int dx = -2; dx < 3; dx++) {
                for (int dy = -2; dy < 3; dy++) {
                    auto x = inds[0] + dx;
                    auto y = inds[1] + dy;
                    if (x < 0 || x >= GRID_LENGTH || y < 0 || y >= GRID_LENGTH) {
                        continue;
                    }
                    GRID[x][y].force += p.mass / p.density * (p.fe * p.fp).determinant() * sigma(p.fe, p.fp) *
                    grad_weighting_func(x, y, p.position[0], p.position[1]);
                }
            }
        }
        //4
        for (int x = 0; x < GRID_LENGTH; x++) {
            for (int y = 0; y < GRID_LENGTH; y++) {
                if (GRID[x][y].mass == 0) {
                    continue;
                }
                auto change_in_velocity = GRID[x][y].force * (1.0 / GRID[x][y].mass) * DELTA_T;
                GRID[x][y].velocity += change_in_velocity;
                // Gravity
                GRID[x][y].velocity[1] -= 9.8 * DELTA_T;
            }
        }
        //5
        //6
        //7
        for (int i = 0; i < N_PARTICLES; i++) {
            auto p = PARTICLES[i];
            auto inds = to_grid(p.position);
            auto grad_v_p = Matrix2d::Zero();
            for (int dx = -2; dx < 3; dx++) {
                for (int dy = -2; dy < 3; dy++) {
                    auto x = inds[0] + dx;
                    auto y = inds[1] + dy;
                    if (x < 0 || x >= GRID_LENGTH || y < 0 || y >= GRID_LENGTH) {
                        continue;
                    }
                    grad_v_p += GRID[x][y].velocity * grad_weighting_func(x, y, p.position[0], p.position[1]).transpose();
                }
            }

            auto fe_hat = (Matrix2d::Identity() + DELTA_T * grad_v_p) * p.fe;
            auto f_n_plus_1 = fe_hat * p.fp;
            auto fe_svd = fe_hat.jacobiSvd(Eigen::ComputeFullU | Eigen::ComputeFullV);
            auto clamped_sing = fe_svd.singularValues().cwiseMin(1 + CRITICAL_STRETCH).cwiseMax(1 - CRITICAL_STRETCH);
            auto clamped_sing_mat = Matrix2d{{clamped_sing[0], 0}, {0, clamped_sing[1]}};

            p.fe = fe_svd.matrixU() * clamped_sing_mat * fe_svd.matrixV().transpose();
            p.fp = fe_svd.matrixV().transpose() * clamped_sing_mat.inverse() * fe_svd.matrixU().transpose() * f_n_plus_1;
        }
        //8
        const double ALPHA = 0.95;
        for (int i = 0; i < N_PARTICLES; i++) {
            auto p = PARTICLES[i];
            auto inds = to_grid(PARTICLES[i].position);
            
            auto v_pic = Vector2d::Zero();
            auto v_flip = p.velocity;
            for (int dx = -2; dx < 3; dx++) {
                for (int dy = -2; dy < 3; dy++) {
                    auto x = inds[0] + dx;
                    auto y = inds[1] + dy;
                    if (x < 0 || x >= GRID_LENGTH || y < 0 || y >= GRID_LENGTH) {
                        continue;
                    }
                    v_pic += (GRID[x][y].velocity * particle_to_grid(Vector2i{x, y}, p.position));
                    v_flip += (GRID[x][y].velocity - GRID[x][y].prev_velocity) * particle_to_grid(Vector2i{x, y}, p.position);
                }
            }
            PARTICLES[i].velocity = (1 - ALPHA) * v_pic + ALPHA * v_flip;
        }
        //9
        //10
        for (int i = 0; i < N_PARTICLES; i++) {
            PARTICLES[i].position += PARTICLES[i].velocity * DELTA_T;
        }
    }
}