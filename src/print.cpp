#include "dataTypes.h"


void printResult(vector<orbitalCurrentResult> & result)
{
    for_each(result.begin(), result.end(), [] (orbitalCurrentResult VCMO)
    {
        printf("Orbital number:    %i   (  %s  ) \n", VCMO.absPosition, VCMO.relPosition.c_str());
        printf("-------------------------------- \n");
        for_each(VCMO.viPoints.begin(), VCMO.viPoints.end(), [] (integrResult point)
        {
            printf("%10.3f %12.2e %12.2e \n", point.voltage, point.current, point.curr_error);
        });
        printf("-------------------------------- \n \n");
    });
}





void printTotalResult(vector<orbitalCurrentResult> &result)
{
    printf(" ************************* Entering Total result calculation module ************************* \n \n");

    //-------------- define and initialize output variables --------------//
    vector <integrResult> totalIntegrResult;

    for_each(result[0].viPoints.begin(), result[0].viPoints.end(), [& totalIntegrResult] (integrResult viPoint)
    {
        totalIntegrResult.push_back(viPoint);
    });
    //--------------------------------------------------------------------//



    //------------------------------- summation over begin+1 to end orbital ------------------------------------//
    for_each(result.begin()+1, result.end(), [& totalIntegrResult] (orbitalCurrentResult VCMO)
    {
        transform(VCMO.viPoints.begin(), VCMO.viPoints.end(), totalIntegrResult.begin(), totalIntegrResult.begin(), [] (integrResult p, integrResult r)
        {
            r.current += p.current;
            r.curr_error += p.curr_error;
            r.trans_prob += p.trans_prob;
            return r;
        });
    });
    //--------------------------------------------------------------------------------------------------------//



    //------------------------------- print Total result ------------------------------------//

    printf("      Total result as summa over all orbitals of orbital window \n");
    printf(" ------------------------------------------------------------------- \n");
    printf(" Voltage (V)   Current (A)   Current error    Transition probability \n");
    printf(" ------------------------------------------------------------------- \n");

    for_each(totalIntegrResult.begin(), totalIntegrResult.end(), [] (integrResult point)
    {
        printf(" %7.3f %15.2e %15.2e %18.2e \n", point.voltage, point.current, point.curr_error, point.trans_prob);
    });
    printf(" ------------------------------------------------------------------- \n \n \n");

    //---------------------------------------------------------------------------------------//




    //------------------------ make and transform vectors to calculate orbital decomposition in total current -----------------------//
    printf(" Orbital decomposition current, current_error and transition probability into total result \n");
    printf(" ----------------------------------------------------------------------------------------- \n");
    printf(" Orbital Number                  Current %%     Current_error %%    transition probability %% \n");
    printf(" ----------------------------------------------------------------------------------------- \n");



    vec tot_currVec(totalIntegrResult.size(), fill::zeros), tot_curr_errorVec(totalIntegrResult.size(), fill::zeros), tot_trans_probVec(totalIntegrResult.size(), fill::zeros);

    for(size_t i = 0; i < totalIntegrResult.size(); i++)
    {
        tot_currVec(i) = totalIntegrResult[i].current;
        tot_curr_errorVec(i) = totalIntegrResult[i].curr_error;
        tot_trans_probVec(i) = totalIntegrResult[i].trans_prob;
    }





    for_each(result.begin(), result.end(), [&tot_currVec, &tot_curr_errorVec, &tot_trans_probVec] (orbitalCurrentResult VCMO)
    {
        vec currentVec(VCMO.viPoints.size(), fill::zeros), curr_errorVec(VCMO.viPoints.size(), fill::zeros), trans_probVec(VCMO.viPoints.size(), fill::zeros);

        //-------------------- fill vectors --------------------//
        for(size_t i = 0; i < VCMO.viPoints.size(); i++)
        {
            currentVec(i) = VCMO.viPoints[i].current;
            curr_errorVec(i) = VCMO.viPoints[i].curr_error;
            trans_probVec(i) = VCMO.viPoints[i].trans_prob;
        }
        //------------------------------------------------------//




        //--------------- normalize all vectors ---------------//
        vec tot_currVec_n = tot_currVec / norm(tot_currVec);
        currentVec = currentVec / norm(tot_currVec);

        vec tot_curr_errorVec_n = tot_curr_errorVec / norm(tot_curr_errorVec);
        curr_errorVec = curr_errorVec / norm(tot_curr_errorVec);

        vec tot_trans_probVec_n = tot_trans_probVec / norm(tot_trans_probVec);
        trans_probVec = trans_probVec / norm(tot_trans_probVec);
        //----------------------------------------------------//

        printf(" %-4u     %-12s ", VCMO.absPosition, VCMO.relPosition.c_str());
        printf("%16.2f  %16.2f  %18.2f \n", dot(currentVec, tot_currVec_n)*100, dot(curr_errorVec, tot_curr_errorVec_n)*100, dot(trans_probVec, tot_trans_probVec_n)*100);


    });

    printf(" ----------------------------------------------------------------------------------------- \n \n");

    printf(" ************************* Out of Total result calculation module ************************** \n \n");
}

