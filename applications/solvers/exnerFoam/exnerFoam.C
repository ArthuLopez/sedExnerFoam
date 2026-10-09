/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

Application
    exnerFoam

Description
    Idealized Exner solver.

\*---------------------------------------------------------------------------*/
//    \heading Solver details
//    The equation is given by:
//
//    \f[
//        \ddt{zb} + \div \left(\vec{qb} \right) = 0
//    \f]

//    Where:
//    \vartable
//        zb         | Bed elevation
//        \vec{qb}   | bedload flux
//    \endvartable

//    \heading Required fields
//    \plaintable
//        zb      | bed elevation [m]
//        qb      | bedload flux [m2/s]
//    \endplaintable


#include "fvCFD.H"
#include "faCFD.H"
#include "dynamicFvMesh.H"
#include "meshTools.H"
#include "pimpleControl.H"
#include "CorrectPhi.H"
#include "localEulerDdtScheme.H"
#include "unitConversion.H"
#include "fvcSmooth.H"
#include "primitivePatchInterpolation.H"
#include "pointMesh.H"
#include "pointPatchField.H"

#include "sedimentBed.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    argList::addNote
    (
        "Scalar transport and incompressible turbulent flow solver."
    );

    #include "addCheckCaseOptions.H"
    #include "setRootCaseLists.H"
    #include "createTime.H"
    #include "createDynamicFvMesh.H"
    #include "createVolFields.H"
    #include "createFaFields.H"

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info<< "\nStarting time loop\n" << endl;

    while (runTime.run())
    {
        ++runTime;

        Info<< "Time = " << runTime.timeName() << nl << endl;

        mesh.controlledUpdate();

        // zb match bed level
        zb = - bed.aMesh().areaCentres() & eg;

        // face normal vectors
        const vectorField& nFaces = bed.aMesh().faceAreaNormals().internalField();

        // compute bedload flux based on H - zb
        forAll(bed.aMesh().areaCentres(), facei)
        {

            // bedload computation
            vector ui = Q[facei] / (H[facei] - zb[facei]);
            scalar magUi = Foam::mag(ui);
            // make ui tangential to bed

	    vector t = ui - (ui & nFaces[facei]) * nFaces[facei];
            t /= Foam::mag(t);

            qsat[facei] = alphaQb * Foam::pow(magUi, betaQb) * t;

	    beta[facei] = Foam::acos(nFaces[facei] & eg);
            // slope direction
            
            if(switchAvalanche)
            {
            	vector slopeDir = eg - nFaces[facei] * (eg & nFaces[facei]);
            	slopeDir /= (Foam::mag(slopeDir) + SMALL);

            	qav[facei] = qav0*slopeDir*Foam::pos(beta[facei] - betaRep)*(Foam::tanh(Foam::tan(beta[facei])) - Foam::tanh(Foam::tan(betaRep)))
                  / (1 - Foam::tanh(Foam::tan(betaRep)));
            }
            else
            {
		qav[facei] = Zero;
            }

        }

        Info << "max Q avalanching: " << Foam::max(Foam::mag(qav)) << endl;
        Info << "max slope: " << radToDeg()*Foam::max(beta) << "°" << endl;
        
	qsat.correctBoundaryConditions();

	if (switchSaturation)
	{

            areaVectorField qsatHat = qsat/(Foam::mag(qsat) + qSmall);

	    areaTensorField outerTensor = qsatHat * qb;

    	    faVectorMatrix saturationEqn
	    (
    		Tsat*fam::ddt(qb)
  		+ Lsat*fac::div(outerTensor)
  		+ fam::Sp(1.0, qb)
 		==
    		qsat
	    );

    	    saturationEqn.solve();
      	}
	else
	{
    	    qb = qsat;
	}

        qb.correctBoundaryConditions();

        phiqb = fac::interpolate(qb) & bed.aMesh().Le();
        phiqav = fac::interpolate(qav) & bed.aMesh().Le();

        // explicit resolution
        faScalarMatrix exnerEqn
        (
            fam::ddt(zb)
            ==
          - fac::div(qb)
          - fac::div(qav)
        );

        exnerEqn.solve();

        dzb = zb - zb.oldTime();

        // for(int edgei=0; edgei < bed.aMesh().nInternalEdges(); edgei++)
        // {
        //     const labelList& owner = bed.aMesh().owner();
        //     Info << "edge " << edgei << "; qb_owner = " << qb[owner[edgei]]
        //         << "; phiqb = " << phiqb[edgei]
        //         << "; dzbi = " << dzb[owner[edgei]] << endl;
        // }

        #include "moveMesh.H"

        if (runTime.writeTime())
        {
            // map areaFields to volFields for vizualisation
            bed.vsm.ref().mapToVolume
            (
                qav,
                qavVf.boundaryFieldRef()
            );
            bed.vsm.ref().mapToVolume
            (
                qb,
                qbVf.boundaryFieldRef()
            );
            runTime.write();

            runTime.printExecutionTime(Info);
        }
    }

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //

