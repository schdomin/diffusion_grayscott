#include "CDomain.h" //ds grid
#include "Timer.h"   //ds time measurement
#include <stdlib.h>  //ds atoi
#include <math.h>    //ds sqrt, etc.
#include <iostream>  //ds cout

int main( int argc, char** argv )
{
    //ds check simple input arguments - CAUTION: the implementation expects real numbers, the simulation will be corrupted if invalid values are entered
    if( 3 != argc )
    {
        //ds inform
        std::cout << "usage: diffusion_pse [Number of grid points 1D] [Final time]" << std::endl;
        return 0;
    }

    //ds start timing
    Timer tmTimer; tmTimer.start( );

    //ds get parameters
    const double dDiffusionCoefficientU( 0.00002 );
    const double dDiffusionCoefficientV( 0.00001 );
    const std::pair< double, double > prBoundaries( -1.0, 1.0 );
    const double dBoundarySize = fabs( prBoundaries.second - prBoundaries.first );
    const unsigned int uNumberOfGridPoints1D( atoi( argv[1] ) );
    const unsigned int uNumberOfParticles( uNumberOfGridPoints1D*uNumberOfGridPoints1D );
    const double dGridPointSpacing( dBoundarySize/( uNumberOfGridPoints1D - 1 ) );
    const double dFinalTime( atoi( argv[2] ) );
    double dTimeStepSize( 1.0 );

    //ds determine timestep size
    if( dDiffusionCoefficientU > dDiffusionCoefficientV )
    {
        dTimeStepSize = 0.5*dGridPointSpacing*dGridPointSpacing/dDiffusionCoefficientU;
    }
    else
    {
        dTimeStepSize = 0.5*dGridPointSpacing*dGridPointSpacing/dDiffusionCoefficientV;
    }

    //ds compute number of steps
    const unsigned int uNumberOfTimeSteps( ceil( dFinalTime/dTimeStepSize ) );

    const double dReactionRateF( 0.03 );
    const double dReactionRateK( 0.062 );

    //ds user information
    std::cout << "\n---------------------------- DIFFUSION PSE GRAY-SCOTT SETUP ----------------------------" << std::endl;
    std::cout << " Diffusion Coefficient u: "  << dDiffusionCoefficientU << std::endl;
    std::cout << " Diffusion Coefficient v: "  << dDiffusionCoefficientV << std::endl;
    std::cout << "           Boundary (2D): [" << prBoundaries.first << ", " << prBoundaries.second << "]" << std::endl;
    std::cout << "Number of Grid Points 1D: "  << uNumberOfGridPoints1D << std::endl;
    std::cout << "     Number of Particles: "  << uNumberOfParticles << std::endl;
    std::cout << "      Grid Point Spacing: "  << dGridPointSpacing << std::endl;
    std::cout << "              Final Time: "  << dFinalTime << std::endl;
    std::cout << "    Number of Time Steps: "  << uNumberOfTimeSteps << std::endl;
    std::cout << "          Time Step Size: "  << dTimeStepSize << std::endl;
    std::cout << "         Reaction rate F: "  << dReactionRateF << std::endl;
    std::cout << "         Reaction rate k: "  << dReactionRateK << std::endl;
    std::cout << "----------------------------------------------------------------------------------------" << std::endl;

    //ds allocate domain (automatically creates initial density distribution)
    Diffusion::CDomain cDomain( dDiffusionCoefficientU, dDiffusionCoefficientV, prBoundaries, dBoundarySize, uNumberOfGridPoints1D, uNumberOfParticles, dGridPointSpacing, dTimeStepSize );

    //ds information
    std::cout << "               Status:  0% done - current time: 0";

    //ds start simulation
    for( unsigned int uCurrentTimeStep = 1; uCurrentTimeStep < uNumberOfTimeSteps+1; ++uCurrentTimeStep )
    {
        //ds calculate percentage done
        const double dPercentageDone( 100.0*uCurrentTimeStep/uNumberOfTimeSteps );

        //ds and time
        const double dCurrentTime( uCurrentTimeStep*dTimeStepSize );

        //ds get a formatted string -> 100% -> 3 digits
        char chBuffer[4];

        //ds fill the buffer
        std::snprintf( chBuffer, 4, "%3.0f", dPercentageDone );

        //ds print info
        std::cout << '\xd';
        std::cout << "               Status: " << chBuffer << "% done - current time: " << dCurrentTime;

        //ds update domain
        cDomain.updateHeatDistributionNumerical( dReactionRateF, dReactionRateK );

        //ds streaming
        cDomain.saveHeatGridToStream( );
    }

    //ds save final mesh to png
    cDomain.saveMeshToPNG( dFinalTime );

    //ds save the streams to a file
    cDomain.writeHeatGridToFile( "bin/simulation.txt", uNumberOfTimeSteps );

    //ds stop timing
    const double dDurationSeconds( tmTimer.stop( ) );

    //ds cause an output ostream
    std::cout << std::endl;
    std::cout << "     Computation time: " << dDurationSeconds << std::endl;
    std::cout << "----------------------------------------------------------------------------------------\n" << std::endl;

    return 0;
}
