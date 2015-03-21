// The FOAM Project // File: util/Ucomponents/Ucomponents.C
/*
-------------------------------------------------------------------------------
 =========         | Application
 \\      /         |
  \\    /          | Name:   Ucomponents
   \\  /           | Family: Utility
    \\/            |
    F ield         | FOAM version: 1.9.5
    O peration     |
    A and          | Copyright (C) 1991-1999 Henry G Weller
    M anipulation  |          All Rights Reserved.
-------------------------------------------------------------------------------
APPLICATION
    

DESCRIPTION
    

AUTHOR
    

-------------------------------------------------------------------------------
*/

#include "fvCFD.H"
#include "Random.H"
#include "mathematicalConstants.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
#   include "addLatestTimeOption.H"
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
    if (Uheader.headerOk())
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
        const scalar betaPlus = 2.0*mathematicalConstant::pi*(1.0/200.0);
        //const scalar sigma = 0.00055;
        const scalar sigma = 0.0002;
        //streamwise wave number: spacing x+ = 500
        const scalar alphaPlus = 2.0*mathematicalConstant::pi*(1.0/800.0);
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
