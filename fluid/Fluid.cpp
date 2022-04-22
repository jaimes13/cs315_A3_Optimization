
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
void Fluid::Create(float w, float h) 
{
	width = w;
	height = h;
	grid_w = (int)(width / FluidSmoothLen);
	grid_h = (int)(height / FluidSmoothLen);

	delete[] gridoffsets;
	gridoffsets = new FluidGridOffset[ grid_w * grid_h ];

	planes.push_back(D3DXVECTOR3(1, 0, 0));
	planes.push_back(D3DXVECTOR3(0, 1, 0));
	planes.push_back(D3DXVECTOR3(-1, 0, width));
	planes.push_back(D3DXVECTOR3(0, -1, height));
}

// Fill a region in the lower left with evenly spaced particles
void Fluid::Fill(float size) 
{
	Clear();

	int w = (int)(size / FluidInitialSpacing);

	// Allocate
	for (int i = 0; i < w*w; i++)
	{
		Particle* p = new Particle;
		particles.push_back(p);
		reconstruction_particles.push_back(nullptr);
	}

	// Populate
	for ( int x = 0 ; x < w ; x++ )
	{
		for ( int y = 0 ; y < w ; y++ )	 
		{
			particles[y*w + x]->pos = D3DXVECTOR2(x * FluidInitialSpacing, Height() - y * FluidInitialSpacing);
			particles[y*w + x]->vel = D3DXVECTOR2(0, 0);
			particles[y*w + x]->acc = D3DXVECTOR2(0, 0);
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
	reconstruction_particles.clear();
}

// Expand the Neighbors list if necessary
// This function is rarely called
__forceinline void Fluid::ExpandNeighbors() 
{
	// Increase the size of the neighbors array because it is full
	neighbors_capacity *= 2;
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
	float cell_count = grid_w * grid_h;
	for( int offset = 0; offset < cell_count; offset++ ) 
	{
		gridoffsets[offset].count = 0;
	}

	// Count the number of particles in each cell
	for( unsigned int particle = 0; particle < particles.size(); particle++ ) 
	{
		// Find where this particle is in the grid
		int p_gx = min(max((int)(particles[particle]->pos.x / FluidSmoothLen), 0), grid_w - 1);
		int p_gy = min(max((int)(particles[particle]->pos.y / FluidSmoothLen), 0), grid_h - 1);
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
		int p_gx = min(max((int)(particles[particle]->pos.x / FluidSmoothLen), 0), grid_w - 1);
		int p_gy = min(max((int)(particles[particle]->pos.y / FluidSmoothLen), 0), grid_h - 1);
		int cell = p_gy * grid_w + p_gx ;
		reconstruction_particles[ gridoffsets[ cell ].offset + gridoffsets[ cell ].count ] = particles[particle];
		gridoffsets[ cell ].count++;
	}
	particles = reconstruction_particles;
}

// Simulation Update
// Build a list of neighbors (particles from adjacent grid locations) for every particle
void Fluid::GetNeighbors() 
{
	// Search radius is the smoothing length
	float h2 = FluidSmoothLen*FluidSmoothLen;

	num_neighbors = 0;
	
	for( unsigned int P = 0; P < particles.size(); P++ )
	{
		// Find where this particle is in the grid
		int p_gx = min(max((int)(particles[P]->pos.x / FluidSmoothLen), 0), grid_w - 1);
		int p_gy = min(max((int)(particles[P]->pos.y / FluidSmoothLen), 0), grid_h - 1);
		int cell = p_gy * grid_w + p_gx ;
		D3DXVECTOR2 pos_P = particles[P]->pos;

		// For every adjacent grid cell (9 cells total for 2D)
		float d_gx_max = ((p_gx < grid_w - 1) ? 1 : 0);
		float d_gy_max = ((p_gy < grid_h - 1) ? 1 : 0);

		for (int d_gx = ((p_gx<1)?0:-1); d_gx <= d_gx_max; d_gx++) 
		{
			for (int d_gy = ((p_gy<1)?0:-1); d_gy <= d_gy_max; d_gy++) 
			{
				// Neighboring cell
				int n_cell = cell + d_gy * grid_w + d_gx; 

				// Loop over ever particle in the neighboring cell
				unsigned int start = gridoffsets[n_cell].offset;
				unsigned int count = gridoffsets[n_cell].count;

				for (int i = 0; i < count; i++)
				{
					unsigned int N = start + i;
					// Only record particle "pairs" once
					if (P > N) 
					{
						// Distance squared
						D3DXVECTOR2 d = pos_P - particles[N]->pos;
						float distsq = d.x * d.x + d.y * d.y;

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
	float constant = (FluidSmoothLen * FluidSmoothLen) * (FluidSmoothLen * FluidSmoothLen) * (FluidSmoothLen * FluidSmoothLen) * FluidWaterMass;

	for( unsigned int particle = 0; particle < particles.size(); particle++ )
	{
		// This is r = 0
		particles[particle]->density = constant;
	}

	// foreach neighboring pair of particles
	for( unsigned int i = 0; i < num_neighbors ; i++ ) 
	{		
		// distance squared
		float r2 = neighbors[i].distsq;
		
		// Density is based on proximity and mass
		// Density is:
		// M_n * W(h, r)
		// Where the smoothing kernel is:
		// The the "Poly6" kernel
		float h2_r2 = FluidSmoothLen * FluidSmoothLen - r2;
		float dens = h2_r2*h2_r2*h2_r2;

		float mass = FluidWaterMass * dens;
		 
		particles[neighbors[i].p]->density += mass;
		particles[neighbors[i].n]->density += mass;
	}

	// Approximate pressure as an ideal compressible gas
	// based on a spring eqation relating the rest density
	for( unsigned int particle = 0 ; particle < particles.size(); ++particle )
	{
		particles[particle]->density *= poly6_coef;
		particles[particle]->pressure = FluidStiff * max(pow(particles[particle]->density / FluidRestDensity, 3) - 1, 0);
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
	float constant = FluidViscosity * lap_vis_coef;

	// foreach neighboring pair of particles
	for( unsigned int i = 0; i < num_neighbors; i++ ) 
	{				
		// Compute force due to pressure and viscosity
		float h_r = FluidSmoothLen - neighbors[i].distsq;
		D3DXVECTOR2 diff = particles[neighbors[i].n]->pos - particles[neighbors[i].p]->pos;

		// Forces is dependant upon the average pressure and the inverse distance
		// Force due to pressure is:
		// 1/rho_p * 1/rho_n * Pavg * W(h, r)
		// Where the smoothing kernel is:
		// The gradient of the "Spikey" kernel
		D3DXVECTOR2 force = (0.5f * (particles[neighbors[i].p]->pressure + particles[neighbors[i].n]->pressure)* grad_spiky_coef * h_r / neighbors[i].distsq ) * diff;
		
		// Viscosity is based on relative velocity
		// Viscosity is:
		// 1/rho_p * 1/rho_n * Vrel * mu * W(h, r)
		// Where the smoothing kernel is:
		// The laplacian of the "Viscosity" kernel
		force += (  constant * (particles[neighbors[i].n]->vel - particles[neighbors[i].p]->vel) );
		
		// Throw in the common (h-r) * 1/rho_p * 1/rho_n
		force *= h_r / (particles[neighbors[i].p]->density * particles[neighbors[i].n]->density);
		
		// Apply force - equal and opposite to both particles
		D3DXVECTOR2 fluid_force = FluidWaterMass * force;
		particles[neighbors[i].p]->acc += fluid_force;
		particles[neighbors[i].n]->acc -= fluid_force;
	}
}

// Simulation Update
// Integration
void Fluid::Integrate( float dt ) 
{
	// Walls
	D3DXVECTOR2 gravity = D3DXVECTOR2(0, 1);
	for( unsigned int particle = 0 ; particle < particles.size() ; ++particle ) 
	{
		// Walls
		size_t size = planes.size();
		for(unsigned int it = 0; it < size; it++)
		{
			float dist = particles[particle]->pos.x * planes[it].x + particles[particle]->pos.y * planes[it].y + planes[it].z;
			particles[particle]->acc += min(dist, 0) * -FluidStaticStiff * D3DXVECTOR2(planes[it].x, planes[it].y );
		}

		// Acceleration
		particles[particle]->acc += gravity;

		// Integration - Euler-Cromer		
		particles[particle]->vel += dt * particles[particle]->acc;
		particles[particle]->pos += dt * particles[particle]->vel;
		particles[particle]->acc = D3DXVECTOR2(0, 0);
	}
}

// Simulation Update
void Fluid::Update( float dt ) 
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
