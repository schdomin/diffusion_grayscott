#ifndef CDOMAIN_H_
#define CDOMAIN_H_

#include <fstream> //ds streaming to file for visualization



namespace Diffusion
{

class CDomain
{

//ds ctor/dtor
public:

    CDomain( const double& p_dDiffusionCoefficientU,
             const double& p_dDiffusionCoefficientV,
             const std::pair< double, double >& p_prBoundaries,
             const double& p_dBoundarySize,
             const unsigned int& p_uNumberOfGridPoints1D,
             const unsigned int& p_uNumberOfParticles,
             const double& p_dGridPointSpacing,
             const double& p_dTimeStepSize );

    ~CDomain( );

//ds attributes
private:

    //ds heat structures for u and v
    double** m_gridHeatU;
    double** m_gridHeatV;

    //ds coefficients
    const double m_dDiffusionCoefficientU;
    const double m_dDiffusionCoefficientV;

    //domain properties
    const std::pair< double, double > m_prBoundaries;
    const double m_dBoundarySize;
    const double m_dGridPointSpacing;
    const unsigned int m_uNumberOfGridPoints1D;
    const unsigned int m_uNumberOfParticles;
    const double m_dTimeStepSize;

    //ds 2d-kernel
    const double m_dEpsilon;
    const double m_dVolP;
    const double m_PSEFactorU;
    const double m_PSEFactorV;

    //ds stream for offline data - needed for the movie and graphs
    std::string m_strLogHeatDistribution;
    std::string m_strLogNorms;

    //ds PNG
    unsigned int m_uNumberOfImagesSaved;

//ds accessors
public:

    void updateHeatDistributionNumerical( const double& p_dReactionRateF, const double& p_dReactionRateK );
    void saveHeatGridToStream( );
    void saveMeshToPNG( const unsigned int& p_uCurrentTimeStep, const unsigned int& p_uRate );
    void writeHeatGridToFile( const std::string& p_strFilename, const unsigned int& p_uNumberOfTimeSteps ) const;

//ds helpers
private:

    //ds returns grid coordinate
    double _onGrid( const unsigned int& p_uIndex ) const;

    void setInitialHeatDistribution( );

    //kernels
    double getKernelEta( const double p_dXp[2], const double p_dXk[2] ) const;
    double getKernelW( const double& p_dLambda ) const;
    double clamp(const double &x, const double &a, const double &b) const;
    unsigned char* getP2M( const unsigned int& p_uMeshSize ) const;

    //ds gray-scott
    double getChiA( const double& p_dX, const double& p_dY ) const;
    double _getNormallyDistributedNumber( ) const;

}; //class CDomain

} //namespace Diffusion



#endif //CDOMAIN_H_
