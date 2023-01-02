/*---------------------------------------------------------------------------*\
Date: August 24, 2006
Author: Eugene de Villiers
Source: http://www.cfd-online.com/Forums/openfoam-solving/58905-les-turbulent-pipe-flow.html#post192079
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
    perturbUCylinder

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
#include "mathematicalConstants.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
#   include "setRootCase.H"

#   include "createTime.H"
#   include "createMesh.H"

    const vectorField centers(mesh.C());

    Info<< "Time = " << runTime.value() << endl;


    IOobject Uheader
    (
        "U",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ
    );

    // Check U exists
    if (Uheader.typeHeaderOk<volVectorField>())
    {
        Info<< "    Reading U" << endl;
        volVectorField U(Uheader, mesh);
        const scalar Retau = 300;
        const scalar d = 7.54/(2000.0);

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

        const scalar utau = Retau*nu.value()/d;
        //wall normal circulation
        const scalar duplus = Ubar.value()[0]*0.5/utau;
        //spanwise wavenumber: spacing z+ = 200
        const scalar betaPlus = 2.0*constant::mathematical::pi*(1.0/200.0);
        //const scalar sigma = 0.00055;
        const scalar sigma = 0.0002;
        //streamwise wave number: spacing x+ = 500
        const scalar alphaPlus = 2.0*constant::mathematical::pi*(1.0/800.0);
        const scalar epsilon = Ubar.value()[0]/20.0;

        forAll(centers, celli)
        {
            scalar& Ux(U[celli].x());
            vector cCenter = centers[celli];
            scalar r = ::sqrt(::sqr(cCenter.y()) + ::sqr(cCenter.z()));
            Ux = 2*mag(Ubar.value())*(1-::sqr(r/d));
            scalar y = d - r;
            r = r*Retau/d;
            y = y*Retau/d;

            scalar theta = ::atan(cCenter.y()/cCenter.z());
            scalar x = cCenter.x()*Retau/d;


            Ux = Ux + (utau*duplus/2.0)
                 *::cos(betaPlus*theta*r) *(y/30)
                 *::exp(-sigma*::sqr(y) + 0.5);
            scalar utheta = epsilon*::sin(alphaPlus*x)*y
                            *::exp(-sigma*::sqr(y));

            vector tangential
            (
                0, cCenter.y(), cCenter.z()
            );
            tangential = tangential ^ vector(1,0,0);
            tangential = tangential/mag(tangential);

            U[celli] = U[celli] + utheta*tangential;

        }

        U.write();
    }
    else
    {
        Info<< "    No U" << endl;
    }

    Info<< endl;


    return(0);
}


// ************************************************************************* //
