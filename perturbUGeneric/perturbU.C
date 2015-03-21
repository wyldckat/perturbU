/*---------------------------------------------------------------------------*\
Date: March 2, 2007
Author: Eugene de Villiers
Source: http://www.cfd-online.com/Forums/openfoam-solving/58905-les-turbulent-pipe-flow-2.html#post192103
-------------------------------------------------------------------------------
License
    This file is a derivative work of OpenFOAM.

    OpenFOAM is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 2 of the License, or (at your
    option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM; if not, write to the Free Software Foundation,
    Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA

Application
    perturbU

Description
    initialise channel velocity with superimposed streamwise streaks.
    To be used to force channelOodles to transition and reach a fully
    developed flow sooner.

    Reads in perturbUDict.

    EdV from paper:
        Schoppa, W. and Hussain, F.
        "Coherent structure dynamics in near wall turbulence",
        Fluid Dynamics Research, Vol 26, pp119-139, 2000.

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "Random.H"
#include "wallDist.H"
#include "wallDistReflection.H"
#include "mathematicalConstants.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
#   include "setRootCase.H"

#   include "createTime.H"
#   include "createMesh.H"

    Info << "perturbU: generalised velocity perturbation implementation for "
         << "the initialisation of ducted well-resolved LES flows." << endl;
    
    //Xdir -> Ubar - streamwise
    //Ydir -> wallReflection vectors
    //Zdir -> cross product Xdir^Ydir
    
    //Ubar and Retau should be set in transportProperties
    //The choice of Retau is not critical as long as it is 
    //approximately the right order of magnitude
    
    //A laminar background velocity profile is assumed
    //with maximum U at h = max(wall distance)
    //A laminar initial profile is essential since wall normal motion
    //of a mean turbulent profile without resolved turbulence will
    //diffuse the perturbations, preventing transition in some cases
    
    
    wallDist yw(mesh);
    const scalar h = max(yw.internalField());
    
    //local yDir
    wallDistReflection reflexVec(mesh);
    
    const volVectorField yDir = reflexVec.n();


    IOobject Uheader
    (
        "U",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ
    );
    Info << "Reading U" << endl;
    volVectorField U(Uheader, mesh);

    IOdictionary transportProperties
    (
        IOobject
        (
            "transportProperties",
            runTime.constant(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    );

    dimensionedScalar nu
    (
        transportProperties.lookup("nu")
    );
    dimensionedVector Ubar
    (
        transportProperties.lookup("Ubar")
    );
    dimensionedScalar Retau
    (
        transportProperties.lookup("Retau")
    );
    
    vector xDir = Ubar.value()/mag(Ubar.value());


    Info << "Re(tau) = " << Retau << endl;
    const scalar utau = Retau.value()*nu.value()/h;
    Info << " u(tau) = " << utau << endl;


    //wall normal circulation
    const scalar duplus = Ubar.value().x()*0.25/utau;
    //spanwise wavenumber: spacing z+ = 200
    const scalar betaPlus = 2.0*mathematicalConstant::pi*(1.0/200.0);
    const scalar sigma = 0.00055;
    //streamwise wave number: spacing x+ = 500
    const scalar alphaPlus = 2.0*mathematicalConstant::pi*(1.0/500.0);
    const scalar epsilon = Ubar.value().x()/200.0;

    const vectorField& centres(mesh.C());

    Random perturbation(1234567);

    forAll(centres, celli)
    {
        //add a small random component to enhance symmetry breaking
        scalar deviation=1.0 + 0.1*perturbation.GaussNormal();
        const vector& cCentre = centres[celli];

        vector zDir = xDir^yDir[celli];
        zDir /= mag(zDir);

        scalar zplus = (cCentre & zDir)*Retau.value()/h;
        scalar yplus = yw[celli]*Retau.value()/h;
        scalar xplus = (cCentre & xDir)*Retau.value()/h;

        // laminar parabolic profile
        U[celli] = 3.0*Ubar.value() * (yw[celli]/h - 0.5*sqr(yw[celli]/h));

        // streak streamwise velocity
        U[celli] +=
            xDir*(utau * duplus/2.0) * (yplus/40.0)
            * Foam::exp(-sigma * Foam::sqr(yplus) + 0.5)
            * Foam::cos(betaPlus*zplus)*deviation;

        // streak spanwise perturbation
        U[celli] +=
            zDir * epsilon
          * Foam::sin(alphaPlus*xplus)
          * yplus
          * Foam::exp(-sigma*Foam::sqr(yplus))
          * deviation;

    }

    U.write();

    Info<< endl;

    return(0);
}


// ************************************************************************* //
