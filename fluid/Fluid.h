
#pragma once

#pragma warning(push, 0)
#include <d3dx9math.h>
#pragma warning(pop)
#include <vector>


// Fluid magic numbers
const float FluidTimestep = 0.005f;
const float FluidSmoothLen = 0.012f;
const float FluidStaticStiff = 3000.0f;
const float FluidRestDensity = 1000.0f;
const float FluidWaterMass = 0.0002f;
const float FluidViscosity = 0.1f;
const float FluidStiff = 200.0f;
const float FluidInitialSpacing = 0.0045f;

/*****************************************************************************/

struct FluidNeighborRecord 
{
	unsigned int p; // Particle Index
	unsigned int n; // Neighbor Index		
	float distsq; // Distance Squared and Non-Squared (After Density Calc)
};

struct FluidGridOffset 
{
	unsigned int offset; // offset into gridindices
	unsigned int count; // number of particles in cell
};

/*****************************************************************************/

struct Particle
{
	D3DXVECTOR2 pos;
	D3DXVECTOR2 vel;
	D3DXVECTOR2 acc;
	float density = 0;
	float pressure = 0;
};

class Fluid 
{
	public:
		/* Common Interface */
		Fluid();
		~Fluid();

		void Create(float w, float h);
		void Fill(float size);
		void Clear();
		void Update(float dt);

		/* Common Data */
		std::vector<Particle*> particles;
		std::vector<Particle*> reconstruction_particles;


		FluidGridOffset * gridoffsets;
		unsigned int neighbors_capacity;
		unsigned int num_neighbors;
		FluidNeighborRecord * neighbors;
		
		unsigned int Size()					{ return particles.size(); }
		unsigned int Step()					{ return step; }
		void Pause( bool p )				{ paused = p; }
		void PauseOnStep( unsigned int p )	{ pause_step = p; }
		float Width()						{ return width; }
		float Height()						{ return height; }

	private:
		
		/* Simulation */
		void UpdateGrid();
		__forceinline void ExpandNeighbors();
		void GetNeighbors();
		void ComputeDensity();
		void SqrtDist();
		void ComputeForce();
		void Integrate(float dt);

	private:
		/* Run State */
		unsigned int step;
		bool paused;
		unsigned int pause_step;

		/* World Size */
		float width;
		float height;
		int grid_w;
		int grid_h;

		/* Coefficients for kernel */
		float poly6_coef;
		float grad_spiky_coef;
		float lap_vis_coef;
		std::vector<D3DXVECTOR3> planes;
};
