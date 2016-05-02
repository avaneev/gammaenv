#include <stdio.h>
#include "CDSPModGammaEnv.h"

int main()
{
	vox :: CDSPModGammaEnv env;
	const double SampleRate = 44100;
	const double Time = 0.050;
	const int SteepCount = 10;
	const int z = (int) floor( Time * 1.5 * SampleRate );
	double tbl[ SteepCount ][ z ];
	int i;
	int j;

	for( j = 0; j < SteepCount; j++ )
	{
		env.Attack = Time;
		env.Release = Time;
		env.AttackDelay = 0.25 * j / ( SteepCount - 1 );
		env.ReleaseDelay = 0.25 * j / ( SteepCount - 1 );
		env.IsInverse = false;
		env.init( SampleRate );
		env.clear( 0.25 );

		for( i = 0; i < z; i++ )
		{
			tbl[ j ][ i ] = env.process( 1.0 );
		}
	}

	for( i = 0; i < z; i += 1 )
	{
		for( j = 0; j < SteepCount - 1; j++ )
		{
			printf( "%f\t", tbl[ j ][ i ]);
		}

		printf( "%f\n", tbl[ j ][ i ]);
	}
}
