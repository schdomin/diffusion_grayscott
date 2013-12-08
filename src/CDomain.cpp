#include "CDomain.h"
#include <math.h>    //ds fabs, etc.
#include <stdlib.h>  //ds rand48
#include "writepng.h"
#include <iostream>



//ds speed
static const double M_PI_SQUARED = M_PI*M_PI;

namespace Diffusion
{

//ds ctor/dtor
CDomain::CDomain( const double& p_dDiffusionCoefficientU,
                  const double& p_dDiffusionCoefficientV,
                  const std::pair< double, double >& p_prBoundaries,
                  const double& p_dBoundarySize,
                  const unsigned int& p_uNumberOfGridPoints1D,
                  const unsigned int& p_uNumberOfParticles,
                  const double& p_dGridPointSpacing,
                  const double& p_dTimeStepSize ) : m_dDiffusionCoefficientU( p_dDiffusionCoefficientU ),
                                                    m_dDiffusionCoefficientV( p_dDiffusionCoefficientV ),
                                                    m_prBoundaries( p_prBoundaries ),
                                                    m_dBoundarySize( p_dBoundarySize ),
                                                    m_uNumberOfGridPoints1D( p_uNumberOfGridPoints1D ),
                                                    m_uNumberOfParticles( p_uNumberOfParticles ),
                                                    m_dGridPointSpacing( p_dGridPointSpacing ),
                                                    m_dTimeStepSize( p_dTimeStepSize ),
                                                    m_dEpsilon( 2*p_dGridPointSpacing ),
                                                    m_dVolP( p_dGridPointSpacing*p_dGridPointSpacing ),
                                                    m_PSEFactorU( p_dDiffusionCoefficientU/( m_dEpsilon*m_dEpsilon )*m_dVolP ),
                                                    m_PSEFactorV( p_dDiffusionCoefficientV/( m_dEpsilon*m_dEpsilon )*m_dVolP ),
                                                    m_strLogHeatDistribution( "" ),
                                                    m_strLogNorms( "" )

{
    //ds allocate memory for the data structures
    m_gridHeatU = new double*[m_uNumberOfGridPoints1D];
    m_gridHeatV = new double*[m_uNumberOfGridPoints1D];

    //ds for each element
    for( unsigned int x = 0; x < m_uNumberOfGridPoints1D; ++x )
    {
        m_gridHeatU[x] = new double[m_uNumberOfGridPoints1D];
        m_gridHeatV[x] = new double[m_uNumberOfGridPoints1D];
    }

    //ds initialize grid
    setInitialHeatDistribution( );
};

CDomain::~CDomain( )
{
    //ds deallocate heat structure
    for( unsigned int x = 0; x < m_uNumberOfGridPoints1D; ++x )
    {
        delete[] m_gridHeatU[x];
        delete[] m_gridHeatV[x];
    }

    delete[] m_gridHeatU;
    delete[] m_gridHeatV;

};

//ds accessors
void CDomain::updateHeatDistributionNumerical( const double& p_dReactionRateF, const double& p_dReactionRateK )
{
    //ds heat change for current time step
    double gridHeatUChangePSE[m_uNumberOfGridPoints1D][m_uNumberOfGridPoints1D];
    double gridHeatVChangePSE[m_uNumberOfGridPoints1D][m_uNumberOfGridPoints1D];

    //ds for all grid points
    for( unsigned int uk = 0; uk < m_uNumberOfGridPoints1D; ++uk )
    {
        for( unsigned int vk = 0; vk < m_uNumberOfGridPoints1D; ++vk )
        {
            //ds get 2d vector of current coordinates
            const double dXk[2] = { _onGrid( uk ), _onGrid( vk ) };

            //ds inner sum of formula
            double dInnerSumU( 0.0 );
            double dInnerSumV( 0.0 );

            //ds loop 20x20
            for( int i = -10; i <= 10; ++i )
            {
                for( int j = -10; j <= 10; ++j )
                {
                    //ds get current indexes up, vp
                    int up = uk + i;
                    int vp = vk + j;

                    //ds if we are not ourself (unsigned overflow no problem here)
                    if( uk != static_cast< unsigned int >( up ) && vk != static_cast< unsigned int >( vp ) )
                    {
                        //ds offset value (required for spaced positions)
                        double dOffsetU = 0.0;
                        double dOffsetV = 0.0;

                        //ds check boundary
                        if( 0 > up                        ){ up += m_uNumberOfGridPoints1D; dOffsetU = -m_dBoundarySize; } //ds moving to negative boundary
                   else if( m_uNumberOfGridPoints1D <= up ){ up -= m_uNumberOfGridPoints1D; dOffsetU = m_dBoundarySize;  } //ds moving to positive boundary
                        if( 0 > vp                        ){ vp += m_uNumberOfGridPoints1D; dOffsetV = -m_dBoundarySize; } //ds moving to negative boundary
                   else if( m_uNumberOfGridPoints1D <= vp ){ vp -= m_uNumberOfGridPoints1D; dOffsetV = m_dBoundarySize;  } //ds moving to positive boundary

                        //ds get 2d vector of current coordinates
                        const double dXp[2] = { ( _onGrid( up ) + dOffsetU ), ( _onGrid( vp ) + dOffsetV ) };

                        //ds get kernel value
                        const double dEta( getKernelEta( dXp, dXk ) );

                        //ds compute inner sum
                        dInnerSumU += ( m_gridHeatU[up][vp] - m_gridHeatU[uk][vk] )*dEta;
                        dInnerSumV += ( m_gridHeatV[up][vp] - m_gridHeatV[uk][vk] )*dEta;
                    }
                }
            }

            //ds heat product
            const double dHeatProduct( m_gridHeatU[uk][vk]*m_gridHeatV[uk][vk]*m_gridHeatV[uk][vk] );

            //ds add final part of formula and save in temporary grid
            gridHeatUChangePSE[uk][vk] = m_PSEFactorU*dInnerSumU - dHeatProduct + p_dReactionRateF*( 1.0 - m_gridHeatU[uk][vk] );
            gridHeatVChangePSE[uk][vk] = m_PSEFactorV*dInnerSumV + dHeatProduct - ( p_dReactionRateF + p_dReactionRateK )*m_gridHeatV[uk][vk];
        }
    }

    //ds copy all computed values to original grid
    for( unsigned int u = 0; u < m_uNumberOfGridPoints1D; ++u )
    {
        for( unsigned int v = 0; v < m_uNumberOfGridPoints1D; ++v )
        {
            m_gridHeatU[u][v] += m_dTimeStepSize*gridHeatUChangePSE[u][v];
            m_gridHeatV[u][v] += m_dTimeStepSize*gridHeatVChangePSE[u][v];
        }
    }
}

void CDomain::saveHeatGridToStream( )
{
    //ds buffer for snprintf
    char chBuffer[16];

    //ds add each element
    for( unsigned int u = 0; u < m_uNumberOfGridPoints1D; ++u )
    {
        for( unsigned int v = 0; v < m_uNumberOfGridPoints1D; ++v )
        {
            //ds get the integrals stream
            std::snprintf( chBuffer, 16, "%f", m_gridHeatU[u][v] );

            //ds add buffer and space to string
            m_strLogHeatDistribution += chBuffer;
            m_strLogHeatDistribution += " ";
        }

        //ds add new line for new row
        m_strLogHeatDistribution += '\n';
    }
}

void CDomain::saveMeshToPNG( const double& p_dCurrentTime )
{
    //ds get mesh size
    const unsigned int uMeshSize( 2*m_uNumberOfGridPoints1D );

    //ds get the data
    unsigned char* chMesh( getP2M( uMeshSize ) );

    //ds construct picture name - buffer for snprintf
    char chBuffer[64];

    //ds format: image_number.png
    std::snprintf( chBuffer, 64, "bin/image_%06.0f.png", p_dCurrentTime );

    //ds create png
    writePNG( chBuffer, uMeshSize, uMeshSize, chMesh );

    //ds free data
    delete chMesh;
}

void CDomain::writeHeatGridToFile( const std::string& p_strFilename, const unsigned int& p_uNumberOfTimeSteps ) const
{
    //ds ofstream object
    std::ofstream ofsFile;

    //ds open the file for writing
    ofsFile.open( p_strFilename.c_str( ), std::ofstream::out );

    //ds if it worked
    if( ofsFile.is_open( ) )
    {
        //ds first dump setup information number of points and time steps
        ofsFile << m_uNumberOfGridPoints1D << " " << p_uNumberOfTimeSteps << "\n" << m_strLogHeatDistribution;
    }

    //ds close the file
    ofsFile.close( );
}

//ds helpers
double CDomain::_onGrid( const unsigned int& p_uIndex ) const
{
    return ( p_uIndex*m_dGridPointSpacing + m_prBoundaries.first );
}

void CDomain::setInitialHeatDistribution( )
{
    //ds loop over all indexi
    for( unsigned int x = 0; x < m_uNumberOfGridPoints1D; ++x )
    {
        for( unsigned int y = 0; y < m_uNumberOfGridPoints1D; ++y )
        {
            //ds get grid coordinates
            const double dX( _onGrid( x ) );
            const double dY( _onGrid( y ) );

            //ds set initial heat values
            m_gridHeatU[x][y] = ( 1.0 - getChiA( dX, dY ) ) + getChiA( dX, dY )*( 1.0/2.0 + _getNormallyDistributedNumber( )/100.0 );
            m_gridHeatV[x][y] = getChiA( dX, dY )*( 1.0/4.0 + _getNormallyDistributedNumber( )/100.0 );
        }
    }
}

double CDomain::getKernelEta( const double p_dXp[2], const double p_dXk[2] ) const
{
    //ds compute distance (using pow for readability)
    const double dDistance = sqrt( ( p_dXp[0] - p_dXk[0] )*( p_dXp[0] - p_dXk[0] ) + ( p_dXp[1] - p_dXk[1] )*( p_dXp[1] - p_dXk[1] ) );

    //ds if we are out of the cutoff
    if( 5*m_dEpsilon < dDistance )
    {
        return 0.0;
    }
    else
    {
        //return kernel function
        return 16.0/( M_PI_SQUARED*( pow( dDistance, 8 ) + 1.0 ) );
    }
}

double CDomain::getKernelW( const double& p_dLambda ) const
{
    //ds get absolute lambda
    const double dLambdaAbs = fabs( p_dLambda );

    //ds cases
    if( 0 <= dLambdaAbs && dLambdaAbs < 1 )
    {
        return ( 1.0 - 5.0/2.0*dLambdaAbs*dLambdaAbs + 3.0/2.0*dLambdaAbs*dLambdaAbs*dLambdaAbs );
    }
    else if( 1 <= dLambdaAbs && dLambdaAbs < 2 )
    {
        return ( 1.0/2.0*( 2.0 - dLambdaAbs )*( 2.0 - dLambdaAbs )*( 1.0 - dLambdaAbs ) );
    }
    else
    {
        return 0.0;
    }
}

double CDomain::clamp(const double &x, const double &a, const double &b) const
{
    return x < a ? a : x > b ? b : x;
}

unsigned char* CDomain::getP2M( const unsigned int& p_uMeshSize ) const
{
    //ds h value
    const double dH( 1.0/p_uMeshSize );

    //ds allocate data array
    unsigned char* chMesh = new unsigned char[p_uMeshSize*p_uMeshSize*4];

    //ds set data - for each data element
    for( unsigned int i = 0; i < p_uMeshSize; ++i )
    {
        for( unsigned int j = 0; j < p_uMeshSize; ++j )
        {
            //ds index in mesh
            unsigned int uIndexMesh = 4.0*( i*p_uMeshSize + j );

            //ds concentrations on mesh
            double dConcentrationU( 0.0 );
            double dConcentrationV( 0.0 );

            //ds for all particles
            for( unsigned int u = 0; u < m_uNumberOfGridPoints1D; ++u )
            {
                for( unsigned int v = 0; v < m_uNumberOfGridPoints1D; ++v )
                {
                    //ds get lambdas
                    const double dLambdai = _onGrid( u )/dH - i;
                    const double dLambdaj = _onGrid( v )/dH - j;

                    //ds get kernel
                    const double dW( getKernelW( dLambdai )*getKernelW( dLambdaj ) );

                    //ds compute concentration
                    dConcentrationU += m_gridHeatU[u][v]*dW;
                    dConcentrationV += m_gridHeatV[u][v]*dW;
                }
            }

            //ds compute difference
            const double dDifference( dConcentrationU-dConcentrationV );

            //ds set the mesh values - SNIPPET
            chMesh[uIndexMesh+0] = clamp( dDifference, 0.0, 1.0 )*0xFFu;  // red
            chMesh[uIndexMesh+1] = clamp( 1.0, 0.0, 1.0 )*0xFFu;          // green
            chMesh[uIndexMesh+2] = clamp( -dDifference, 0.0, 1.0 )*0xFFu; // blue
            chMesh[uIndexMesh+3] = 0xFFu;                                 // alpha (opacity, 0xFFu = 255)
        }
    }

    return chMesh;
}

double CDomain::getChiA( const double& p_dX, const double& p_dY ) const
{
    //ds case
    if( (-0.2 <= p_dX && p_dX <= 0.2 ) && (-0.2 <= p_dY && p_dY <= 0.2 ) )
    {
        return 1.0;
    }
    else
    {
        return 0.0;
    }

}

double CDomain::_getNormallyDistributedNumber( ) const
{
    //ds calculate the uniform number first [0,1]
    const double dUniformNumber( drand48( ) );

    //ds return the normal one
    return sqrt( -2*log( dUniformNumber ) )*cos( 2*M_PI*dUniformNumber );
}

} //namespace Diffusion
