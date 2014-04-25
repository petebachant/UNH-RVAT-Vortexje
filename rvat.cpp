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

class Blade : public LiftingSurface
{
public:
    // Constructor:
    Blade()
    {
        // Create blade:
        LiftingSurfaceBuilder surface_builder(*this);
        
        const int n_points_per_airfoil = 32;
        const int n_airfoils = 10; // Was 21
        
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
        const int n_layers = 10; // was 21
        
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

class Walls : public Body
{
public:
    Walls(string   id) : Body(id)
    {
        double extrude_length = 4.0;
        
        // Create floor
        Vector3d start_point(-2, -1.83, -1.22);
        int span_direction = 1; // Span in the y direction (a floor)
        double span_length = 3.66;
        
        Surface *floor = new Wall(start_point,
                                  span_direction,
                                  span_length,
                                  extrude_length);
        add_non_lifting_surface(*floor);
        allocated_surfaces.push_back(floor);
        
        // Create top
        start_point(2) = 1.22;
        span_direction = 1; // Span in the y direction (a floor)
        span_length = 3.66;
        
        Surface *top = new Wall(start_point,
                                span_direction,
                                span_length,
                                extrude_length);
        add_non_lifting_surface(*top);
        allocated_surfaces.push_back(top);
        
        // Create right wall (looking downstream)
        start_point(2) = -1.22;
        span_direction = 2; // Span in the z direction
        span_length = 2.44;
        
        Surface *rightwall = new Wall(start_point,
                                      span_direction,
                                      span_length,
                                      extrude_length);
        add_non_lifting_surface(*rightwall);
        allocated_surfaces.push_back(rightwall);
        
        // Create left wall (looking downstream)
        start_point(1) = 1.83;
        span_direction = 2; // Span in the z direction
        span_length = 2.44;
        
        Surface *leftwall = new Wall(start_point,
                                     span_direction,
                                     span_length,
                                     extrude_length);
        add_non_lifting_surface(*leftwall);
        allocated_surfaces.push_back(leftwall);
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
              M_PI / 6.0,
              TIP_SPEED_RATIO * WIND_VELOCITY / MILL_RADIUS);
              
    Walls walls(string("walls"));
    
    // Set up solver:
    Solver solver("rvat-log");
    solver.add_body(vawt);
    solver.add_body(walls);
    
    Vector3d freestream_velocity(WIND_VELOCITY, 0, 0);
    solver.set_freestream_velocity(freestream_velocity);
    
    double fluid_density = 1000.0;
    solver.set_fluid_density(fluid_density);
    
    // Set up file format for logging:
    VTKSurfaceWriter surface_writer;
    
    // Set up velocity field writer
    VTKFieldWriter field_writer;
    double dx = 0.1;
    double dy = 0.05;
    double dz = 0.125;
    double x_min = 1.0; 
    double y_min = -1.5; 
    double z_min = 0;
    double x_max = 1.1; 
    double y_max = 1.5; 
    double z_max = 0.625;
    mkdir("rvat-log/velocity", S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
    
    // Log shaft moments:
    ofstream f;
    f.open("rvat-log/performance.txt");
    
    // Run simulation:
    double t = 0.0;
    double dt = 0.02; // Was 0.0033
    int step_number = 0;
    
    solver.initialize_wakes(dt);
    while (t < 10) {
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
        velfilename << "rvat-log/velocity/" << t << ".vtk";
        field_writer.write_velocity_field(solver, velfilename.str(), 
                                          x_min, x_max,
                                          y_min, y_max,
                                          z_min, z_max,
                                          dx, dy, dz);
        
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
