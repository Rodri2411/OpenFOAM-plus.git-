/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2017 OpenCFD Ltd.
     \\/     M anipulation  |
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

\*---------------------------------------------------------------------------*/

#include "objectRegistry.H"
#include "Time.H"

template<class Type>
bool Foam::functionObjects::fieldAverageItem::calculateMeanField
(
    const objectRegistry& obr
) const
{
    if (!mean_)
    {
        return false;
    }

    const Type* baseFieldPtr = obr.lookupObjectPtr<Type>(fieldName_);

    if (!baseFieldPtr)
    {
        return false;
    }

    const Type& baseField = *baseFieldPtr;
    Type& meanField = obr.lookupObjectRef<Type>(meanFieldName_);

    switch (windowType_)
    {
        case windowType::NONE:
        {
            scalar dt = this->dt(obr.time().deltaTValue());
            scalar Dt = this->Dt();
            scalar beta = dt/Dt;

            meanField = (1 - beta)*meanField + beta*baseField;

            break;
        }
        case windowType::APPROXIMATE:
        {
            scalar dt = this->dt(obr.time().deltaTValue());
            scalar Dt = this->Dt();
            scalar beta = dt/Dt;

            if (Dt - dt >= window_)
            {
                beta = dt/window_;
            }

            meanField = (1 - beta)*meanField + beta*baseField;

            break;
        }
        case windowType::EXACT:
        {
            switch (base_)
            {
                case baseType::ITER:
                {
                    // Uniform time step - can use simplified algorithm
                    // Note: stores an additional old time field, but only
                    //       needs to do 1 field lookup

                    label n = windowTimes_.size();

                    if (n <= round(window_))
                    {
                        scalar beta = 1.0/scalar(n);
                        meanField = (1 - beta)*meanField + beta*baseField;
                    }
                    else
                    {
                        const Type& lastField =
                            obr.lookupObject<Type>(windowFieldNames_.first());
                        meanField += (baseField - lastField)/scalar(n - 1);
                    }

                    break;
                }
                case baseType::TIME:
                {
                    if (windowTimes_.size() < 2)
                    {
                        meanField = 1*baseField;
                    }
                    else
                    {
                        meanField = 0*baseField;
                        FIFOStack<scalar>::const_iterator timeIter =
                            windowTimes_.begin();
                        FIFOStack<word>::const_iterator nameIter =
                            windowFieldNames_.begin();

                        const Type* w0Ptr =
                            obr.lookupObjectPtr<Type>(nameIter());

                        ++nameIter;
                        scalar t0 = timeIter();
                        ++timeIter;

                        for
                        (
                            ;
                            timeIter != windowTimes_.end();
                            ++timeIter, ++nameIter
                        )
                        {
                            const Type* wPtr =
                                obr.lookupObjectPtr<Type>(nameIter());
                            const scalar t = timeIter();
                            const Type& w = *wPtr;
                            const Type& w0 = *w0Ptr;

                            if (t0 > window_)
                            {
                                scalar tStar = max(0, window_ - t);
                                meanField +=
                                    ((w0 - w)/(t0 - t)*tStar + 2*w)*tStar;
                            }
                            else
                            {
                                meanField += (w + w0)*(t0 - t);
                            }

                            w0Ptr = wPtr;
                            t0 = t;
                        }

                        meanField /= 2*min(windowTimes_.first(), window_);
                    }

                    break;
                }
                default:
                {
                    FatalErrorInFunction
                        << "Unhandled baseType enumeration "
                        << baseTypeNames_[base_]
                        << abort(FatalError);
                }
            }

            break;
        }
        default:
        {
            FatalErrorInFunction
                << "Unhandled windowType enumeration "
                << windowTypeNames_[windowType_]
                << abort(FatalError);
        }
    }

    return true;
}


template<class Type1, class Type2>
bool Foam::functionObjects::fieldAverageItem::calculatePrime2MeanField
(
    const objectRegistry& obr
) const
{
    if (!prime2Mean_)
    {
        return false;
    }

    const Type1* baseFieldPtr = obr.lookupObjectPtr<Type1>(fieldName_);

    if (!baseFieldPtr)
    {
        return false;
    }

    const Type1& baseField = *baseFieldPtr;
    const Type1& meanField = obr.lookupObject<Type1>(meanFieldName_);

    Type2& prime2MeanField =
        obr.lookupObjectRef<Type2>(prime2MeanFieldName_);

    switch (windowType_)
    {
        case windowType::NONE:
        {
            scalar dt = this->dt(obr.time().deltaTValue());
            scalar Dt = this->Dt();
            scalar beta = dt/Dt;

            prime2MeanField =
                (1 - beta)*prime2MeanField
              + beta*sqr(baseField)
              - sqr(meanField);

            break;
        }
        case windowType::APPROXIMATE:
        {
            scalar dt = this->dt(obr.time().deltaTValue());
            scalar Dt = this->Dt();
            scalar beta = dt/Dt;

            if (Dt - dt >= window_)
            {
                beta = dt/window_;
            }

            prime2MeanField =
                (1 - beta)*prime2MeanField
              + beta*sqr(baseField)
              - sqr(meanField);

            break;
        }
        case windowType::EXACT:
        {
            FIFOStack<scalar>::const_iterator timeIter =
                windowTimes_.begin();
            FIFOStack<word>::const_iterator nameIter =
                windowFieldNames_.begin();

            if (windowFieldNames_.size() < 2)
            {
                prime2MeanField = sqr(baseField - meanField);
            }
            else
            {
                switch (base_)
                {
                    case baseType::ITER:
                    {
                        prime2MeanField = 0*prime2MeanField;
                        label n = min(windowFieldNames_.size(), window_);
                        DebugVar(windowFieldNames_.size());
                        DebugVar(window_);
                        if (windowFieldNames_.size() > window_)
                        {
                            ++nameIter;
                            ++timeIter;
                        }

                        for
                        (
                            ;
                            nameIter != windowFieldNames_.end();
                            ++nameIter
                        )
                        {
                            const Type1* wPtr =
                                obr.lookupObjectPtr<Type1>(nameIter());
                            const Type1& w = *wPtr;

                            prime2MeanField += sqr(w - meanField);
                        }

                        prime2MeanField /= n;

                        break;
                    }
                    case baseType::TIME:
                    {
                        prime2MeanField = 0*prime2MeanField;
                        const Type1* w0Ptr =
                            obr.lookupObjectPtr<Type1>(nameIter());

                        ++nameIter;
                        scalar t0 = timeIter();
                        ++timeIter;
                        for
                        (
                            ;
                            timeIter != windowTimes_.end();
                            ++timeIter, ++nameIter
                        )
                        {
                            const Type1* w1Ptr =
                                obr.lookupObjectPtr<Type1>(nameIter());
                            const scalar t1 = timeIter();

                            const Type1& w0 = *w0Ptr;
                            const Type1& w1 = *w1Ptr;

                            if (t0 > window_)
                            {
                                scalar tStar = max(0, window_ - t1);
                                prime2MeanField +=
                                    tStar
                                   *sqr
                                    (
                                        0.5*((w1 - w0)/(t1 - t0)*tStar + 2*w1)
                                      - meanField
                                    );
                            }
                            else
                            {
                                prime2MeanField +=
                                    (t0 - t1)*(sqr(0.5*(w0 + w1) - meanField));
                            }

                            w0Ptr = w1Ptr;
                            t0 = t1;
                        }

                        prime2MeanField /= min(windowTimes_.first(), window_);

                        break;
                    }
                    default:
                    {
                        FatalErrorInFunction
                            << "Unhandled baseType enumeration "
                            << baseTypeNames_[base_]
                            << abort(FatalError);
                    }
                }
            }

            break;
        }
        default:
        {
            FatalErrorInFunction
                << "Unhandled windowType enumeration "
                << windowTypeNames_[windowType_]
                << abort(FatalError);
        }
    }

    return true;
}


// ************************************************************************* //
