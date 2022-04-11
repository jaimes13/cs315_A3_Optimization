
#include "Fluid.h"

// Zero the fluid simulation member variables for sanity
Fluid::Fluid() 
{
	step = 0;
	paused = false;
	pause_step = 0xFFFFFFFF;

	width = 0;
	height = 0;
	grid_w = 0;
	grid_h = 0;

	gridindices = NULL;
	gridoffsets = NULL;
	num_neighbors = 0;
	// If this value is too small, ExpandNeighbors will fix it
	neighbors_capacity = 263 * 1200;
	neighbors = new FluidNeighborRecord[ neighbors_capacity ];

	// Precompute kernel coefficients
	// See "Particle-Based Fluid Simulation for Interactive Applications"
	// "Poly 6" Kernel - Used for Density
	poly6_coef = 315.0f / (64.0f * D3DX_PI * pow(FluidSmoothLen, 9));
	// Gradient of the "Spikey" Kernel - Used for Pressure
	grad_spiky_coef = -45.0f / (D3DX_PI * pow(FluidSmoothLen, 6));
	// Laplacian of the "Viscosity" Kernel - Used for Viscosity
	lap_vis_coef = 45.0f / (D3DX_PI * pow(FluidSmoothLen, 6));
}

// Destructor
Fluid::~Fluid() 
{
	Clear();
	delete[] gridoffsets; gridoffsets = NULL;
	num_neighbors = 0;
	neighbors_capacity = 0;
	delete[] neighbors; neighbors = neighbors;
}

// Create the fluid simulation
// width/height is the simulation world maximum size
void Fluid::Create(double w, double h) 
{
	width = w;
	height = h;
	grid_w = (int)(width / FluidSmoothLen);
	grid_h = (int)(height / FluidSmoothLen);

	delete[] gridoffsets;
	gridoffsets = new FluidGridOffset[ grid_w * grid_h ];
}

Particle* Fluid::particle_at(std::size_t index)
{
	if (index > particles.size())
		throw "out of bounds";
	
	auto it = particles.begin();
	for (; index != 0; it++, index--);

	return *it;
}

// Fill a region in the lower left with evenly spaced particles
void Fluid::Fill(double size) 
{
	Clear();

	int w = (int)(size / FluidInitialSpacing);

	// Allocate
	gridindices = new unsigned int[ w * w ];
	for (int i = 0; i < w*w; i++)
	{
		Particle* p = new Particle;
		particles.push_back(p);
	}

	// Populate
	for ( int x = 0 ; x < w ; x++ )
	{
		for ( int y = 0 ; y < w ; y++ )	 
		{
			particle_at(y*w + x)->pos = D3DXVECTOR2(x * FluidInitialSpacing, Height() - y * FluidInitialSpacing);
			particle_at(y*w + x)->vel = D3DXVECTOR2(0, 0);
			particle_at(y*w + x)->acc = D3DXVECTOR2(0, 0);
			gridindices[ y*w+x ] = 0;
		}
	}
}

// Remove all particles
void Fluid::Clear() 
{
	step = 0;

	for (auto p : particles)
		delete p;
	particles.clear();

	delete[] gridindices; gridindices = NULL;
}

// Expand the Neighbors list if necessary
// This function is rarely called
__forceinline void Fluid::ExpandNeighbors() 
{
	// Increase the size of the neighbors array because it is full
	neighbors_capacity += 20;
	FluidNeighborRecord* new_neighbors = new FluidNeighborRecord[ neighbors_capacity ];
	memcpy( new_neighbors, neighbors, sizeof(FluidNeighborRecord) * num_neighbors );
	delete[] neighbors;
	neighbors = new_neighbors;
}

// Simulation Update
// Build the grid of neighbors
// Imagine an evenly space grid.  All of our neighbors will be
// in our cell and the 8 adjacent cells
void Fluid::UpdateGrid() 
{
	// Cell size is the smoothing length

	// Clear the offsets
	for( int offset = 0; offset < (grid_w * grid_h); offset++ ) 
	{
		gridoffsets[offset].count = 0;
	}

	// Count the number of particles in each cell
	for( unsigned int particle = 0; particle < particles.size(); particle++ ) 
	{
		// Find where this particle is in the grid
		int p_gx = min(max((int)(particle_at(particle)->pos.x * (1.0 / FluidSmoothLen)), 0), grid_w - 1);
		int p_gy = min(max((int)(particle_at(particle)->pos.y * (1.0 / FluidSmoothLen)), 0), grid_h - 1);
		int cell = p_gy * grid_w + p_gx ;
		gridoffsets[ cell ].count++;
	}

	// Prefix sum all of the cells
	unsigned int sum = 0;
	for( int offset = 0; offset < (grid_w * grid_h); offset++ ) 
	{
		gridoffsets[offset].offset = sum;
		sum += gridoffsets[offset].count;
		gridoffsets[offset].count = 0;
	}

	// Insert the particles into the grid
	for( unsigned int particle = 0; particle < particles.size(); particle++ ) 
	{
		// Find where this particle is in the grid
		int p_gx = min(max((int)(particle_at(particle)->pos.x * (1.0 / FluidSmoothLen)), 0), grid_w - 1);
		int p_gy = min(max((int)(particle_at(particle)->pos.y * (1.0 / FluidSmoothLen)), 0), grid_h - 1);
		int cell = p_gy * grid_w + p_gx ;
		gridindices[ gridoffsets[ cell ].offset + gridoffsets[ cell ].count ] = particle;
		gridoffsets[ cell ].count++;
	}
}

// Simulation Update
// Build a list of neighbors (particles from adjacent grid locations) for every particle
void Fluid::GetNeighbors() 
{
	// Search radius is the smoothing length
	double h2 = FluidSmoothLen*FluidSmoothLen;

	num_neighbors = 0;
	
	for( unsigned int P = 0; P < particles.size(); P++ )
	{
		// Find where this particle is in the grid
		int p_gx = min(max((int)(particle_at(P)->pos.x * (1.0f / FluidSmoothLen)), 0), grid_w - 1);
		int p_gy = min(max((int)(particle_at(P)->pos.y * (1.0f / FluidSmoothLen)), 0), grid_h - 1);
		int cell = p_gy * grid_w + p_gx ;
		D3DXVECTOR2 pos_P = particle_at(P)->pos;

		// For every adjacent grid cell (9 cells total for 2D)
		for (int d_gy = ((p_gy<1)?0:-1); d_gy <= ((p_gy<grid_h-1)?1:0); d_gy++) 
		{
			for (int d_gx = ((p_gx<1)?0:-1); d_gx <= ((p_gx<grid_w-1)?1:0); d_gx++) 
			{
				// Neighboring cell
				int n_cell = cell + d_gy * grid_w + d_gx; 

				// Loop over ever particle in the neighboring cell
				unsigned int* start = gridindices + gridoffsets[n_cell].offset;
				unsigned int* end = start + gridoffsets[n_cell].count;

				for ( ; start != end ; ++start) 
				{
					unsigned int N = *start;
					// Only record particle "pairs" once
					if (P > N) 
					{
						// Distance squared
						D3DXVECTOR2 d = pos_P - particle_at(N)->pos;
						double distsq = d.x * d.x + d.y * d.y;

						// Check that the particle is within the smoothing length
						if (distsq < h2) 
						{
							if ( num_neighbors >= neighbors_capacity ) 
							{
								ExpandNeighbors();
							}
							// Record the ID of the two particles
							// And record the squared distance
							FluidNeighborRecord& record = neighbors[ num_neighbors ];
							record.p = P;
							record.n = N;
							record.distsq = distsq;
							num_neighbors++;
						}
					}
				}
			}
		}
	}
}

// Simulation Update
// Compute the density for each particle based on its neighbors within the smoothing length
void Fluid::ComputeDensity() 
{
	for( unsigned int particle = 0; particle < particles.size(); particle++ )
	{
		// This is r = 0
		particle_at(particle)->density = (FluidSmoothLen * FluidSmoothLen) * (FluidSmoothLen * FluidSmoothLen) * (FluidSmoothLen * FluidSmoothLen) * FluidWaterMass;
	}

	// foreach neighboring pair of particles
	for( unsigned int i = 0; i < num_neighbors ; i++ ) 
	{		
		// distance squared
		double r2 = neighbors[i].distsq;
		
		// Density is based on proximity and mass
		// Density is:
		// M_n * W(h, r)
		// Where the smoothing kernel is:
		// The the "Poly6" kernel
		double h2_r2 = FluidSmoothLen * FluidSmoothLen - r2;
		double dens = h2_r2*h2_r2*h2_r2;

		double P_mass = FluidWaterMass;
		double N_mass = FluidWaterMass;
		 
		particle_at(neighbors[i].p)->density += N_mass * dens;
		particle_at(neighbors[i].n)->density += P_mass * dens;
	}

	// Approximate pressure as an ideal compressible gas
	// based on a spring eqation relating the rest density
	for( unsigned int particle = 0 ; particle < particles.size(); ++particle )
	{
		particle_at(particle)->density *= poly6_coef;
		particle_at(particle)->pressure = FluidStiff * max(pow(particle_at(particle)->density / FluidRestDensity, 3) - 1, 0);
	}
}

// Simulation Update
// Perform a batch of sqrts to turn distance squared into distance
void Fluid::SqrtDist() 
{
	for( unsigned int i = 0; i < num_neighbors; i++ ) 
	{
		neighbors[i].distsq = sqrt(neighbors[i].distsq);
	}
}

// Simulation Update
// Compute the forces based on the Navier-Stokes equations for laminer fluid flow
// Follows is lots more voodoo
void Fluid::ComputeForce() 
{
	// foreach neighboring pair of particles
	for( unsigned int i = 0; i < num_neighbors; i++ ) 
	{				
		// Compute force due to pressure and viscosity
		double h_r = FluidSmoothLen - neighbors[i].distsq;
		D3DXVECTOR2 diff = particle_at(neighbors[i].n)->pos - particle_at(neighbors[i].p)->pos;

		// Forces is dependant upon the average pressure and the inverse distance
		// Force due to pressure is:
		// 1/rho_p * 1/rho_n * Pavg * W(h, r)
		// Where the smoothing kernel is:
		// The gradient of the "Spikey" kernel
		D3DXVECTOR2 force = (0.5f * (particle_at(neighbors[i].p)->pressure + particle_at(neighbors[i].n)->pressure)* grad_spiky_coef * h_r / neighbors[i].distsq ) * diff;
		
		// Viscosity is based on relative velocity
		// Viscosity is:
		// 1/rho_p * 1/rho_n * Vrel * mu * W(h, r)
		// Where the smoothing kernel is:
		// The laplacian of the "Viscosity" kernel
		force += ( (FluidViscosity * lap_vis_coef) * (particle_at(neighbors[i].n)->vel - particle_at(neighbors[i].p)->vel) );
		
		// Throw in the common (h-r) * 1/rho_p * 1/rho_n
		force *= h_r * 1.0f / (particle_at(neighbors[i].p)->density * particle_at(neighbors[i].n)->density);
		
		// Apply force - equal and opposite to both particles
		particle_at(neighbors[i].p)->acc += FluidWaterMass * force;
		particle_at(neighbors[i].n)->acc -= FluidWaterMass * force;
	}
}

// Simulation Update
// Integration
void Fluid::Integrate( double dt ) 
{
	// Walls
	std::list<D3DXVECTOR3> planes;
	planes.push_back( D3DXVECTOR3(1, 0, 0) );
	planes.push_back( D3DXVECTOR3(0, 1, 0) );
	planes.push_back( D3DXVECTOR3(-1, 0, width) );
	planes.push_back( D3DXVECTOR3(0, -1, height) );

	D3DXVECTOR2 gravity = D3DXVECTOR2(0, 1);
	for( unsigned int particle = 0 ; particle < particles.size() ; ++particle ) 
	{
		// Walls
		for( auto it = planes.begin(); it != planes.end(); it++ )
		{
			double dist = particle_at(particle)->pos.x * (*it).x + particle_at(particle)->pos.y * (*it).y + (*it).z;
			particle_at(particle)->acc += min(dist, 0) * -FluidStaticStiff * D3DXVECTOR2( (*it).x, (*it).y );
		}

		// Acceleration
		particle_at(particle)->acc += gravity;

		// Integration - Euler-Cromer		
		particle_at(particle)->vel += dt * particle_at(particle)->acc;
		particle_at(particle)->pos += dt * particle_at(particle)->vel;
		particle_at(particle)->acc = D3DXVECTOR2(0, 0);
	}
}

// Simulation Update
void Fluid::Update( double dt ) 
{
	// Pause runs the simulation standing still for profiling
	if( paused || step == pause_step ) { dt = 0.0f; }
	else { step++; }

	// Create neighbor information
	UpdateGrid();
	GetNeighbors();

	// Calculate the forces for all of the particles
	ComputeDensity();
	SqrtDist();
	ComputeForce();

	// And integrate
	Integrate(dt);
}
