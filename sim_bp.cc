#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "sim_bp.h"
#include <math.h>
#include <algorithm>

/*  argc holds the number of command line arguments
    argv[] holds the commands themselves

    Example:-
    sim bimodal 6 gcc_trace.txt
    argc = 4
    argv[0] = "sim"
    argv[1] = "bimodal"
    argv[2] = "6"
    ... and so on
*/
int main(int argc, char *argv[])
{
    FILE *FP;               // File handler
    char *trace_file;       // Variable that holds trace file name;
    bp_params params;       // look at sim_bp.h header file for the the definition of struct bp_params
    char outcome;           // Variable holds branch outcome
    unsigned long int addr; // Variable holds the address read from input file

    if (!(argc == 4 || argc == 5 || argc == 7))
    {
        printf("Error: Wrong number of inputs:%d\n", argc - 1);
        exit(EXIT_FAILURE);
    }

    params.bp_name = argv[1];

    // strtoul() converts char* to unsigned long. It is included in <stdlib.h>
    if (strcmp(params.bp_name, "bimodal") == 0) // Bimodal
    {
        if (argc != 4)
        {
            printf("Error: %s wrong number of inputs:%d\n", params.bp_name, argc - 1);
            exit(EXIT_FAILURE);
        }
        params.M2 = strtoul(argv[2], NULL, 10);
        params.M1 = 0;
        params.N = 0;
        params.K = 0;
        trace_file = argv[3];
        printf("COMMAND\n%s %s %lu %s\n", argv[0], params.bp_name, params.M2, trace_file);
    }
    else if (strcmp(params.bp_name, "gshare") == 0) // Gshare
    {
        if (argc != 5)
        {
            printf("Error: %s wrong number of inputs:%d\n", params.bp_name, argc - 1);
            exit(EXIT_FAILURE);
        }
        params.M1 = strtoul(argv[2], NULL, 10);
        params.M2 = 0;
        params.K = 0;
        params.N = strtoul(argv[3], NULL, 10);
        trace_file = argv[4];
        printf("COMMAND\n%s %s %lu %lu %s\n", argv[0], params.bp_name, params.M1, params.N, trace_file);
    }
    else if (strcmp(params.bp_name, "hybrid") == 0) // Hybrid
    {
        if (argc != 7)
        {
            printf("Error: %s wrong number of inputs:%d\n", params.bp_name, argc - 1);
            exit(EXIT_FAILURE);
        }
        params.K = strtoul(argv[2], NULL, 10);
        params.M1 = strtoul(argv[3], NULL, 10);
        params.N = strtoul(argv[4], NULL, 10);
        params.M2 = strtoul(argv[5], NULL, 10);
        trace_file = argv[6];
        printf("COMMAND\n%s %s %lu %lu %lu %lu %s\n", argv[0], params.bp_name, params.K, params.M1, params.N, params.M2, trace_file);
    }
    else
    {
        printf("Error: Wrong branch predictor name:%s\n", params.bp_name);
        exit(EXIT_FAILURE);
    }

    // Open trace_file in read mode
    FP = fopen(trace_file, "r");
    if (FP == NULL)
    {
        // Throw error and exit if fopen() failed
        printf("Error: Unable to open file %s\n", trace_file);
        exit(EXIT_FAILURE);
    }

    char str[2];
    int Numberofindex = 0;
    int Numberofindex2 = 0;
    if (strcmp(params.bp_name, "bimodal") == 0)
    {
        Numberofindex2 = pow(2, params.M2);
    }
    if (strcmp(params.bp_name, "gshare") == 0)
    {
        Numberofindex = pow(2, params.M1);
    }
    if (strcmp(params.bp_name, "hybrid") == 0)
    {
        Numberofindex = pow(2, params.M1);
        Numberofindex2 = pow(2, params.M2);
    }
    Numberofindex = pow(2, params.M1);
    Numberofindex2 = pow(2, params.M2);

    int *array = new int[Numberofindex2];
    for (int i = 0; i < Numberofindex2; i++)
    {
        array[i] = 2;
    }

    int *array1 = new int[Numberofindex];
    for (int i = 0; i < Numberofindex; i++)
    {
        array1[i] = 2;
    }

    int numChooserEntries = pow(2, params.K);
    int *chooserTable = new int[numChooserEntries];
    for (int i = 0; i < numChooserEntries; i++)
    {
        chooserTable[i] = 1; // Initialize to 1
    }

    int numberofpredictions = 0;
    int numberofmispredictions = 0;
    float mispredictionrate = 0;
    int NumberofindexBimodal = pow(2, params.M2);
    int NumberofindexGshare = pow(2, params.M1);
    int GlobalHistoryIndex = pow(2, params.N);
    int history = 0;
    history = history & ((1 << params.N) - 1);

    while (fscanf(FP, "%lx %s", &addr, str) != EOF)
    {
        addr = addr >> 2;

        numberofpredictions++;

        outcome = str[0];
        if (outcome == 't')
        {
            if (strcmp(params.bp_name, "bimodal") == 0)
            {
                unsigned int mask = (1 << params.M2) - 1;
                unsigned int Index = addr & mask;
                int get = array[Index];
                if (get < 2)
                {
                    numberofmispredictions++;
                }
                get = get + 1;
                if (get > 3)
                {
                    get = 3;
                }
                array[Index] = get;
            }
            else if (strcmp(params.bp_name, "gshare") == 0)
            {
                unsigned int mask = (1 << params.M1) - 1;
                unsigned int AddressIndex = addr & mask;
                unsigned int firstNBits = AddressIndex >> (params.M1 - params.N);
                int result = firstNBits ^ history;
                unsigned int remainingBits = addr & ((1 << (params.M1 - params.N)) - 1);
                unsigned int finalAddress = (result << (params.M1 - params.N)) | remainingBits;
                history = history >> 1;
                history = history | (1 << (params.N - 1));
                history = history & ((1 << params.N) - 1);
                int get = array1[finalAddress];
                if (get < 2)
                {
                    numberofmispredictions++;
                }
                get = get + 1;
                if (get > 3)
                {
                    get = 3;
                }
                array1[finalAddress] = get;
            }
            else if (strcmp(params.bp_name, "hybrid") == 0)
            {
                // Step 1: Bimodal prediction
                unsigned int bimodal_mask = (1 << params.M2) - 1;
                unsigned int bimodal_index = addr & bimodal_mask;
                int bimodal_prediction = array[bimodal_index];

                // Step 2: Gshare prediction
                unsigned int gshare_mask = (1 << params.M1) - 1;
                unsigned int address_index = addr & gshare_mask;
                unsigned int first_n_bits = address_index >> (params.M1 - params.N);
                int result = first_n_bits ^ history;
                unsigned int remaining_bits = addr & ((1 << (params.M1 - params.N)) - 1);
                unsigned int final_address = (result << (params.M1 - params.N)) | remaining_bits;
                int gshare_prediction = array1[final_address];

                // Step 3: Chooser table prediction
                unsigned int chooser_mask = (1 << params.K) - 1;
                unsigned int chooser_index = (addr)&chooser_mask;

                int chooser_value = chooserTable[chooser_index];

                // Step 4: Choose which predictor to use
                int final_prediction;
                if (chooser_value >= 2)
                {
                    final_prediction = gshare_prediction; // Use gshare if chooser value is >= 2
                    array1[final_address] = std::min(final_prediction + 1, 3);
                }
                else
                {
                    final_prediction = bimodal_prediction; // Use bimodal if chooser value is < 2
                    array[bimodal_index] = std::min(final_prediction + 1, 3);
                }

                // Step 5: Update misprediction count
                if ((final_prediction < 2))
                {
                    numberofmispredictions++;
                }
                // Step 7: Always update the gshare's global history
                history = history >> 1;
                history = history | (1 << (params.N - 1));
                history = history & ((1 << params.N) - 1);

                // Step 8: Update the chooser table based on the accuracy of the predictors
                bool gshare_correct = (gshare_prediction >= 2);
                bool bimodal_correct = (bimodal_prediction >= 2);

                if (gshare_correct && !bimodal_correct)
                {
                    chooserTable[chooser_index] = std::min(chooserTable[chooser_index] + 1, 3);
                }
                else if (bimodal_correct && !gshare_correct)
                {
                    chooserTable[chooser_index] = std::max(chooserTable[chooser_index] - 1, 0);
                }
            }
        }
        else if (outcome == 'n')
        {
            if (strcmp(params.bp_name, "bimodal") == 0)
            {
                unsigned int mask = (1 << params.M2) - 1;
                unsigned int Index = addr & mask;
                int get = array[Index];
                if (get > 1)
                {
                    numberofmispredictions++;
                }
                get = get - 1;
                if (get < 0)
                {
                    get = 0;
                }
                array[Index] = get;
            }
            else if (strcmp(params.bp_name, "gshare") == 0)
            {
                unsigned int mask = (1 << params.M1) - 1;
                unsigned int AddressIndex = addr & mask;
                unsigned int firstNBits = AddressIndex >> (params.M1 - params.N);
                int result = firstNBits ^ history;
                unsigned int remainingBits = addr & ((1 << (params.M1 - params.N)) - 1);
                unsigned int finalAddress = (result << (params.M1 - params.N)) | remainingBits;
                history = history >> 1;
                history = history | (0 << (params.N - 1));
                history = history & ((1 << params.N) - 1);
                int get = array1[finalAddress];
                if (get > 1)
                {
                    numberofmispredictions++;
                }
                get = get - 1;
                if (get < 0)
                {
                    get = 0;
                }
                array1[finalAddress] = get;
            }
            else if (strcmp(params.bp_name, "hybrid") == 0)
            {

                // Step 1: Bimodal prediction
                unsigned int bimodal_mask = (1 << params.M2) - 1;
                unsigned int bimodal_index = addr & bimodal_mask;
                int bimodal_prediction = array[bimodal_index];

                // Step 2: Gshare prediction
                unsigned int gshare_mask = (1 << params.M1) - 1;
                unsigned int address_index = addr & gshare_mask;
                unsigned int first_n_bits = address_index >> (params.M1 - params.N);

                int result = first_n_bits ^ history;
                unsigned int remaining_bits = addr & ((1 << (params.M1 - params.N)) - 1);
                unsigned int final_address = (result << (params.M1 - params.N)) | remaining_bits;
                int gshare_prediction = array1[final_address];

                // Step 3: Chooser table prediction
                unsigned int chooser_mask = (1 << params.K) - 1;
                unsigned int chooser_index = (addr)&chooser_mask;
                int chooser_value = chooserTable[chooser_index];

                // Step 4: Choose which predictor to use
                int final_prediction;
                if (chooser_value >= 2)
                {
                    final_prediction = gshare_prediction; // Use gshare if chooser value is >= 2
                }
                else
                {
                    final_prediction = bimodal_prediction; // Use bimodal if chooser value is < 2
                }

                // Step 5: Update misprediction count
                if ((final_prediction < 2 && outcome == 't') || (final_prediction >= 2 && outcome == 'n'))
                {
                    numberofmispredictions++;
                }
                // Step 8: Update the chooser table based on the accuracy of the predictors
                bool gshare_correct = (gshare_prediction < 2);
                bool bimodal_correct = (bimodal_prediction < 2);

                if (gshare_correct && !bimodal_correct)
                {
                    chooserTable[chooser_index] = std::min(chooserTable[chooser_index] + 1, 3);
                }
                else if (bimodal_correct && !gshare_correct)
                {
                    chooserTable[chooser_index] = std::max(chooserTable[chooser_index] - 1, 0);
                }

                // Step 6: Update the chosen predictor based on the outcome
                if (chooser_value >= 2)
                {
                    // printf("test");
                    gshare_prediction = std::max(gshare_prediction - 1, 0);
                    array1[final_address] = gshare_prediction;
                }
                else
                {
                    // Update bimodal predictor
                    bimodal_prediction = std::max(bimodal_prediction - 1, 0);
                    array[bimodal_index] = bimodal_prediction;
                }

                history = history >> 1;
                history = history | (0 << (params.N - 1));
                history = history & ((1 << params.N) - 1);
            }
        }
    }

    printf("OUTPUT\n");
    printf("number of predictions:    %d\n", numberofpredictions);
    printf("number of mispredictions: %d\n", numberofmispredictions);
    printf("misprediction rate:       %10.2f%%\n", (double(numberofmispredictions) / numberofpredictions) * 100);
    if (strcmp(params.bp_name, "gshare") == 0)
    {
        printf("FINAL GSHARE CONTENTS\n");
        for (int i = 0; i < Numberofindex; i++)
        {
            printf("%-2d %d\n", i, array1[i]);
        }
    }

    else if (strcmp(params.bp_name, "bimodal") == 0)
    {
        printf("FINAL BIMODAL CONTENTS\n");
        for (int i = 0; i < Numberofindex2; i++)
        {
            printf("%-2d %d\n", i, array[i]);
        }
    }

    else if (strcmp(params.bp_name, "hybrid") == 0)
    {
        printf("FINAL CHOOSER CONTENTS\n");
        for (int i = 0; i < numChooserEntries; i++)
        {
            printf("%-2d %d\n", i, chooserTable[i]);
        }

        printf("FINAL GSHARE CONTENTS\n");
        for (int i = 0; i < Numberofindex; i++)
        {
            printf("%-2d %d\n", i, array1[i]);
        }
        printf("FINAL BIMODAL CONTENTS\n");
        for (int i = 0; i < Numberofindex2; i++)
        {
            printf("%-2d %d\n", i, array[i]);
        }
    }

    return 0;
}
