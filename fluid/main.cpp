
#include "Fluid.h"

#pragma warning(push, 0)
#include <windows.h>
#pragma warning(pop)


Fluid * FLUID = NULL;

void InitializeDirectX();
int RunDirectX();

// Frame Update Function - Called Every Frame
void CALLBACK OnFrameMove( double /*fTime*/, float /*fElapsedTime*/, void* /*pUserContext*/ )
{
	FLUID->Update( FluidTimestep );
}

// main
INT WINAPI wWinMain( HINSTANCE, HINSTANCE, LPWSTR, int )
{
	FLUID = new Fluid();
	FLUID->Create( 1.6f, 1.2f );

/*#ifdef _DEBUG
//#error This runs really slow in Debug.  Understand?  Comment out this line, then change the next line to 0.2f
#endif*/
	FLUID->Fill( 0.2f );

	InitializeDirectX();

	int retval = RunDirectX();

	delete FLUID; FLUID = NULL;
	return retval;
}
