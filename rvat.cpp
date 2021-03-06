//
// Adapted from sample by Jorn Baayen <jorn.baayen@baayen-heinz.com>
//

#include <vortexje/solver.hpp>
#include <vortexje/lifting-surface-builder.hpp>
#include <vortexje/shape-generators/airfoils/naca4-airfoil-generator.hpp>
#include <vortexje/shape-generators/ellipse-generator.hpp>
#include <vortexje/surface-writers/vtk-surface-writer.hpp>
#include <vortexje/boundary-layers/dummy-boundary-layer.hpp>
#include <vortexje/empirical-wakes/ramasamy-leishman-wake.hpp>
#include <vortexje/field-writers/vtk-field-writer.hpp>

#include <iostream>
#include <fstream>
#include <sys/stat.h>

using namespace std;
using namespace Eigen;
using namespace Vortexje;

#define N_BLADES        3
#define MILL_RADIUS     0.5
#define TIP_SPEED_RATIO 1.9
#define WIND_VELOCITY   1.0
#define INCLUDE_TOWER
#define INCLUDE_WALLS
#define END_TIME		6.0

#ifdef INCLUDE_WALLS
    string save_dir("rvat-log-walls");
#else
    string save_dir("rvat-log-free");
#endif

class Blade : public LiftingSurface
{
public:
    // Constructor:
    Blade()
    {
        // Create blade:
        LiftingSurfaceBuilder surface_builder(*this);
        
        const int n_points_per_airfoil = 32;
        const int n_airfoils = 16; // Was 21
        
        const double chord = 0.14;
        const double span = 1.0;
        
        int trailing_edge_point_id;
        vector<int> prev_airfoil_nodes;
        
        vector<vector<int> > node_strips;
        vector<vector<int> > panel_strips;
        
        for (int i = 0; i < n_airfoils; i++) {
            vector<Vector3d, Eigen::aligned_allocator<Vector3d> > airfoil_points =
                NACA4AirfoilGenerator::generate(0, 0, 0.20, true, chord, n_points_per_airfoil, trailing_edge_point_id);
            for (int j = 0; j < (int) airfoil_points.size(); j++)
                airfoil_points[j](2) += i * span / (double) (n_airfoils - 1);
                
            vector<int> airfoil_nodes = surface_builder.create_nodes_for_points(airfoil_points);
            node_strips.push_back(airfoil_nodes);
            
            if (i > 0) {
                vector<int> airfoil_panels = surface_builder.create_panels_between_shapes(airfoil_nodes, prev_airfoil_nodes, trailing_edge_point_id);
                panel_strips.push_back(airfoil_panels);
            }
                
            prev_airfoil_nodes = airfoil_nodes;
        }

        surface_builder.finish(node_strips, panel_strips, trailing_edge_point_id);
        
        // Translate and rotate into the canonical coordinate system:
        Vector3d translation(-chord / 2.0, 0.0, -span / 2.0);
        translate(translation);
        
        rotate(Vector3d::UnitZ(), -M_PI / 2.0);
    }
};

class Tower : public Surface
{
public:
    // Constructor:
    Tower()
    {
        // Create cylinder:      
        SurfaceBuilder surface_builder(*this);
        
        const double r = 0.095/2;
        const double h = 1.256;
        
        const int n_points = 32; // was 32
        const int n_layers = 16; // was 21
        
        vector<int> prev_nodes;
        
        for (int i = 0; i < n_layers; i++) {
            vector<Vector3d, Eigen::aligned_allocator<Vector3d> > points =
                EllipseGenerator::generate(r, r, n_points);
            for (int j = 0; j < (int) points.size(); j++)
                points[j](2) += i * h / (double) (n_layers - 1);
                 
            vector<int> nodes = surface_builder.create_nodes_for_points(points);
            
            if (i > 0)
                vector<int> airfoil_panels = surface_builder.create_panels_between_shapes(nodes, prev_nodes);
                
            prev_nodes = nodes;
        }

        surface_builder.finish();
        
        // Translate into the canonical coordinate system:
        Vector3d translation(0.0, 0.0, -h / 2.0);
        translate(translation);
    }
};

// Create a rectangular wall, extruded along x direction
class Wall : public Surface
{
public:
    // Constructor:
    Wall(Vector3d start_point,
         int span_direction, 
         double span_length,
         double extrude_length)
    { 
        SurfaceBuilder surface_builder(*this);
        
        const int n_points = 4;
        const int n_layers = 4;
        
        vector<int> prev_nodes;
        
        // Loop through and move points along extrude direction
        for (int i = 0; i < n_layers; i++) {
            
            // Create points along line
            vector<Vector3d, Eigen::aligned_allocator<Vector3d> > points;
            for (int n = 0; n < n_points; n++){
                Vector3d point = start_point;
                point(span_direction) += n * span_length / (double) (n_points - 1);
                points.push_back(point);
            }
                
            // Move points along extrude direction (x)
            for (int j = 0; j < (int) points.size(); j++)
                points[j](0) += i * extrude_length / (double) (n_layers - 1);
                 
            vector<int> nodes = surface_builder.create_nodes_for_points(points);
            
            if (i > 0)
                vector<int> airfoil_panels = surface_builder.create_panels_between_shapes(nodes, prev_nodes);
                
            prev_nodes = nodes;
        }

        surface_builder.finish();
        
        // Translate into the canonical coordinate system:
        //Vector3d translation(0.0, 0.0, -h / 2.0);
        //translate(translation);
    }
};

// Create a rectangular wall, extruded along x direction
class RectangularTube : public Surface
{
public:
    // Constructor:
    RectangularTube(double height,
                    double width,
                    double extrude_length,
                    bool external_flow = true)
    { 
        SurfaceBuilder surface_builder(*this);
        
        const int n_points = 12;
        const int n_layers = 12;
        
        const double min_x = -3.0;
        
        vector<int> prev_nodes;
        
        // Loop through and move points along extrude direction
        for (int i = 0; i < n_layers; i++) {
            
            // Create points vector
            vector<Vector3d, Eigen::aligned_allocator<Vector3d> > points;
            
            // Lower horizontal line
            for (int n = 0; n < n_points - 1; n++){
                Vector3d point(min_x, -width/2, -height/2);
                point(1) += n * width / (double) (n_points - 1);
                points.push_back(point);
            }
            
            // -y side vertical line
            for (int n = 0; n < n_points - 1; n++){
                Vector3d point(min_x, width/2, -height/2);
                point(2) += n * height / (double) (n_points - 1);
                points.push_back(point);
            }
            
            // Upper horizontal line
            for (int n = 0; n < n_points - 1; n++){
                Vector3d point(min_x, width/2, height/2);
                point(1) -= n * width / (double) (n_points - 1);
                points.push_back(point);
            }
            
            // +y side vertical line
            for (int n = 0; n < n_points - 1; n++){
                Vector3d point(min_x, -width/2, height/2);
                point(2) -= n * height / (double) (n_points - 1);
                points.push_back(point);
            }
                
            // Move points along extrude direction (x)
            for (int j = 0; j < (int) points.size(); j++)
                points[j](0) += i * extrude_length / (double) (n_layers - 1);
                 
            vector<int> nodes = surface_builder.create_nodes_for_points(points);
            
            if (i > 0) {          
                vector<int> airfoil_panels;
                
                if (external_flow)
                    airfoil_panels = surface_builder.create_panels_between_shapes(prev_nodes, nodes);
                else
                    airfoil_panels = surface_builder.create_panels_between_shapes(nodes, prev_nodes);
            }
                
            prev_nodes = nodes;
        }

        surface_builder.finish();
        
    }
};

class Walls : public Body
{
public:
    Walls(string   id) : Body(id)
    {
        double extrude_length = 10.0;
        double height = 2.44;
        double width = 3.66;
        
        Surface *tube = new RectangularTube(height,
                                            width,
					    extrude_length,
                                            false);
        add_non_lifting_surface(*tube);
        allocated_surfaces.push_back(tube);


        Surface *tube_outer = new RectangularTube(height*1.1,
                                                  width*1.1,
						  extrude_length,
                                                  true);
        add_non_lifting_surface(*tube_outer);
        allocated_surfaces.push_back(tube_outer);

    } 
};


class VAWT : public Body
{
public:
    double rotor_radius;
    
    // Constructor:
    VAWT(string   id,
         double   rotor_radius,
         int      n_blades,
         Vector3d position,
         double   theta_0,
         double   dthetadt) :
         Body(id), rotor_radius(rotor_radius)
    {
        // Initialize kinematics:
        this->position = position;
        this->velocity = Vector3d(0, 0, 0);
        this->attitude = AngleAxis<double>(theta_0, Vector3d::UnitZ());
        this->rotational_velocity = Vector3d(0, 0, dthetadt);
        
#ifdef INCLUDE_TOWER
        // Initialize tower:
        Surface *tower = new Tower();
        add_non_lifting_surface(*tower);
        
        allocated_surfaces.push_back(tower);
#endif
        
        // Initialize blades:
        for (int i = 0; i < n_blades; i++) {
            Blade *blade = new Blade();
            
            Vector3d translation(rotor_radius, 0, 0);
            blade->translate(translation);
            
            double theta_blade = theta_0 + 2 * M_PI / n_blades * i;
            blade->rotate(Vector3d::UnitZ(), theta_blade);
            
            blade->translate(position);
            
            BoundaryLayer *boundary_layer = new DummyBoundaryLayer();
            
            Wake *wake = new RamasamyLeishmanWake(*blade);
            
            add_lifting_surface(*blade, *boundary_layer, *wake);
            
            allocated_boundary_layers.push_back(boundary_layer);
            allocated_surfaces.push_back(blade);
            allocated_surfaces.push_back(wake);
        }
    }

    // Rotate:
    void
    rotate(double dt)
    { 
        // Compute new kinematic state:
        Quaterniond new_attitude = AngleAxis<double>(rotational_velocity(2) * dt, Vector3d::UnitZ()) * attitude;
        set_attitude(new_attitude);
    }
};

int
main (int argc, char **argv)
{    
    // Set simulation parameters for dynamic wake convection:
    Parameters::convect_wake = true;
    
    // Set up VAWT:
    Vector3d position(0, 0, 0);
    
    VAWT vawt(string("turbine"),
              MILL_RADIUS,
              N_BLADES,
              position,
              0.0,
              TIP_SPEED_RATIO * WIND_VELOCITY / MILL_RADIUS);
    
    // Set up solver:
    Solver solver(save_dir);
    solver.add_body(vawt);
    
    // Include walls if defined
#ifdef INCLUDE_WALLS
    Walls walls(string("walls"));
    solver.add_body(walls);
#endif
    
    Vector3d freestream_velocity(WIND_VELOCITY, 0, 0);
    solver.set_freestream_velocity(freestream_velocity);
    
    double fluid_density = 1000.0;
    solver.set_fluid_density(fluid_density);
    
    // Set up file format for logging:
    VTKSurfaceWriter surface_writer;
    
    // Set up velocity field writer
    VTKFieldWriter field_writer;
    int nx = 1;
    int ny = 41;
    int nz = 21;
    double x_min = 1.0; 
    double y_min = -1.83; 
    double z_min = -1.22;
    double x_max = 1.2; 
    double y_max = 1.83; 
    double z_max = 1.22;
    string save_subdir = save_dir + string("/velocity");
    mkdir(save_subdir.c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
    
    // Log shaft moments:
    ofstream f;
    save_subdir = save_dir + string("/performance.txt");
    f.open(save_subdir.c_str());
    
    // Run simulation: 
    double t = 0.0; 
    double dt = 0.02; // Was 0.0033
    int step_number = 0;
    
    solver.initialize_wakes(dt);
    while (t < END_TIME) {
        // Solve:
        solver.solve(dt);
        
        // Log coefficients:
        solver.log(step_number, surface_writer);
        
        // Log shaft moment:
		cout << "Computing power coefficient for t = " << t << endl;
        Vector3d M = solver.moment(vawt, position);
		double ct = M(2)/(0.5*fluid_density*1.0*MILL_RADIUS*WIND_VELOCITY*WIND_VELOCITY);
		double cp = ct*TIP_SPEED_RATIO;
		cout << "C_P = " << cp << endl;
        f << t << "    " << cp << endl;
        
        // Write velocity field:
        stringstream velfilename;
        velfilename << save_dir << "/velocity/" << t << ".vtk";
        #ifdef INCLUDE_WALLS
        field_writer.write_velocity_field_surf(solver, velfilename.str(), 
                                               x_min, x_max,
                                               y_min, y_max,
                                               z_min, z_max,
                                               nx, ny, nz,
                                               walls.non_lifting_surfaces[0]->surface);
        #else
        field_writer.write_velocity_field(solver, velfilename.str(), 
                                          x_min, x_max,
                                          y_min, y_max,
                                          z_min, z_max,
                                          nx, ny, nz);
        #endif
        
        // Rotate blades:
        vawt.rotate(dt);
        
        // Update wakes:
        solver.update_wakes(dt);
        
        // Step time:
        t += dt;
        step_number++;
    }
    
    // Close shaft log file:
    f.close();
    
    // Done:
    return 0;
}
