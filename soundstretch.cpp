
#include <stdexcept>
#include <stdio.h>
#include <string.h>
#include <time.h>
#include <math.h>
#include <vector>
#include <iostream>
#include <cfloat>
#include "WavFile.h"
#include "SoundTouch.h"
#include "PeakFinder.h"
#include "dywapitchtrack.h"

using namespace soundtouch;
using namespace std;

typedef struct _tempoidx {
    int idx;
    double tempo;
} tempoidx;

#if _WIN32
    #include <io.h>
    #include <fcntl.h>

    // Macro for Win32 standard input/output stream support: Sets a file stream into binary mode
    #define SET_STREAM_TO_BIN_MODE(f) (_setmode(_fileno(f), _O_BINARY))
#else
    // Not needed for GNU environment... 
    #define SET_STREAM_TO_BIN_MODE(f) {}
#endif

// frame size for interpolation
#define IFRAME_SIZE 20

#define WINDOW_SIZE 31

uint16_t efficient_peak_detector(uint16_t value)
{
    static uint16_t wbuffer[WINDOW_SIZE] = {0U};
    static uint16_t wr_idx = 0U;
    static uint16_t peak_value = 0U;
    static uint16_t peak_value_idx = 0U;

    if (value >= peak_value)
    {
        peak_value = value;            /* New peak value, so record it */
        peak_value_idx = wr_idx;
        wbuffer[wr_idx] = value;
    }
    else
    {
        wbuffer[wr_idx] = value;        /* Not a new peak value, so just store it */
        if (wr_idx == peak_value_idx)    /* Have we over written the peak value ? */
        {
            /*  Yes, so we need to do a brute force search to find the new
                maximum. Note that for efficiency reasons, if there are multiple
                values of the new peak value, then we want to chose the one
                whose index value is as far away as possible from the current index */
            uint16_t idx;
            uint16_t cnt;

            for (cnt = 0U, idx = wr_idx, peak_value = 0U; cnt < WINDOW_SIZE; ++cnt)
            {
                if (wbuffer[idx] >= peak_value)
                {
                    peak_value = wbuffer[idx];    /* Record new peak */
                    peak_value_idx = idx;
                }
                if (++idx >= WINDOW_SIZE)
                {
                    idx = 0;
                }
            }
        }
    }
    if (++wr_idx >= WINDOW_SIZE)
    {
        wr_idx = 0;
    }

    return peak_value;
}

#define STOPPER 0                                      /* Smaller than any datum */
#define    MEDIAN_FILTER_SIZE    (13)

uint16_t median_filter(uint16_t datum)
{
 struct pair
 {
   struct pair   *point;                              /* Pointers forming list linked in sorted order */
   uint16_t  value;                                   /* Values to sort */
 };
 static struct pair buffer[MEDIAN_FILTER_SIZE] = {0}; /* Buffer of nwidth pairs */
 static struct pair *datpoint = buffer;               /* Pointer into circular buffer of data */
 static struct pair small = {NULL, STOPPER};          /* Chain stopper */
 static struct pair big = {&small, 0};                /* Pointer to head (largest) of linked list.*/

 struct pair *successor;                              /* Pointer to successor of replaced data item */
 struct pair *scan;                                   /* Pointer used to scan down the sorted list */
 struct pair *scanold;                                /* Previous value of scan */
 struct pair *median;                                 /* Pointer to median */
 uint16_t i;

 if (datum == STOPPER)
 {
   datum = STOPPER + 1;                             /* No stoppers allowed. */
 }

 if ( (++datpoint - buffer) >= MEDIAN_FILTER_SIZE)
 {
   datpoint = buffer;                               /* Increment and wrap data in pointer.*/
 }

 datpoint->value = datum;                           /* Copy in new datum */
 successor = datpoint->point;                       /* Save pointer to old value's successor */
 median = &big;                                     /* Median initially to first in chain */
 scanold = NULL;                                    /* Scanold initially null. */
 scan = &big;                                       /* Points to pointer to first (largest) datum in chain */

 /* Handle chain-out of first item in chain as special case */
 if (scan->point == datpoint)
 {
   scan->point = successor;
 }
 scanold = scan;                                     /* Save this pointer and   */
 scan = scan->point ;                                /* step down chain */

 /* Loop through the chain, normal loop exit via break. */
 for (i = 0 ; i < MEDIAN_FILTER_SIZE; ++i)
 {
   /* Handle odd-numbered item in chain  */
   if (scan->point == datpoint)
   {
     scan->point = successor;                      /* Chain out the old datum.*/
   }

   if (scan->value < datum)                        /* If datum is larger than scanned value,*/
   {
     datpoint->point = scanold->point;             /* Chain it in here.  */
     scanold->point = datpoint;                    /* Mark it chained in. */
     datum = STOPPER;
   };

   /* Step median pointer down chain after doing odd-numbered element */
   median = median->point;                       /* Step median pointer.  */
   if (scan == &small)
   {
     break;                                      /* Break at end of chain  */
   }
   scanold = scan;                               /* Save this pointer and   */
   scan = scan->point;                           /* step down chain */

   /* Handle even-numbered item in chain.  */
   if (scan->point == datpoint)
   {
     scan->point = successor;
   }

   if (scan->value < datum)
   {
     datpoint->point = scanold->point;
     scanold->point = datpoint;
     datum = STOPPER;
   }

   if (scan == &small)
   {
     break;
   }

   scanold = scan;
   scan = scan->point;
 }
 return median->value;
}

void PitchShift(float pitchShift, long numSampsToProcess, long fftFrameSize, long osamp, float sampleRate, float *indata, float *outdata);

int main(int argc, char* argv[])
{
    long bufferLengthFrames = 512;
    if (argc < 6)
    {
        printf("./meow input.wav 2 2 pattern.wav [pattern0.wav pattern1.wav ... ] output2.wav");
        return -1;
    }

    char *voiceFileName = argv[1];//"input.wav";
    float prettifyAmp = atof( argv[2]); // 2
    float prettifyPitch = atof( argv[3]); // 2

//    char *patternFileName = argv[4]; //"pattern.wav";
    std::vector<std::string> patternFileNames;
    int idx = 4;
    for ( ; idx<argc-1; idx++)
    {
        patternFileNames.push_back(argv[idx]);
//        printf("%s\n", argv[idx]);
    }

    char *outFileName = argv[argc-1];//"output.wav";

    amppitch prevdatas;
    float minAMP = FLT_MAX;
    float maxAMP = 0;
    float summAMP = 0;
    float summPitch = 0;

    float pitchFactor = 0;
    float ampFactor = 0;
    float prevAmpFactor = 0;

    float beginOfWord = 0;
    float endOfWord = 0;
    int beginIdx = 0;
    float stretchFactor = 0;

    int errPitches = 0;

    int insamples;

    std::srand(std::time(0));

    dywapitchtracker voicetracker;
    dywapitch_inittracking(&voicetracker);

    dywapitchtracker patterntracker;
    dywapitch_inittracking(&patterntracker);

    float *voiceData = new float[bufferLengthFrames];

    // Open input file
    WavInFile *infile = new WavInFile(voiceFileName);

    if ((infile->getSampleRate() != 44100) || (infile->getNumChannels() != 1))
    {
        printf("Sorry, this example processes mono audio sampled at 44100Hz for voice.wav.\n");
        return -1;
    }

    for (idx = 0; idx<(int)patternFileNames.size(); idx ++)
    {
        // Open input file
        WavInFile *patternfile = new WavInFile(patternFileNames[idx].c_str());

        if ((patternfile->getSampleRate() != 44100) || (patternfile->getNumChannels() != 1))
        {
           printf("Sorry, this example processes mono audio sampled at 44100Hz.\n");
           delete patternfile;
           return -1;
        }
        delete patternfile;
    }

    uint lenMS = infile->getLengthMS();
    // collect pitches
    std::vector<amppitch> pitches;
    std::vector<int> peaks;

    float prevsum = 0;
    for(;;)
    {
        // Read data from input file
        insamples = infile->read(voiceData, bufferLengthFrames);
        if (infile->eof())
        {
            break;
        }
        amppitch voicedatas = dywapitch_computepitch(&voicetracker, voiceData, 0, bufferLengthFrames);
        pitches.push_back(voicedatas);

        float CUTOFF = 44100/4;
        float RC = 1.0/(CUTOFF*2*3.14);
        float dt = 1.0/44100;
        float alpha = RC/(RC + dt);
        float filteredArray[insamples];

        filteredArray[0] = voiceData[0];
        for (int fidx = 1; fidx<insamples; fidx++){
            filteredArray[fidx] = alpha * (filteredArray[fidx-1] + voiceData[fidx] - voiceData[fidx-1]);
        }

        float summ = 0;
        for (int sidx=0; sidx<insamples; sidx++)
        {
            summ += fabs(filteredArray[sidx]);
        }

        if (((prevsum > 1) && (summ < 1)) || ((prevsum < 1) && (summ > 1)))
        {
//            std::cout << "end of word at " << lenMS/905.0*idx/1000.0 << std::endl;
            peaks.push_back(idx);
        }
        idx++;
        prevsum = summ;
    }

    peaks.push_back(idx);
    delete infile;
    ///////////////////////////////////////////////////// INFILE END

    WavOutFile outFile(outFileName, 44100, 16, 1);

    float deltaMS = lenMS/float(pitches.size());
    std::cout << "pitches stores " << int(pitches.size()) << " numbers. " << lenMS << " " << deltaMS <<"\n";

    insamples = 0;

    int pitchesSize = int(pitches.size());

    // AVG for pitches
    for(idx=0; idx<pitchesSize;idx++)
    {
        amppitch val = pitches[idx];
        if (val.pitch != 0)
        {
            summPitch +=val.pitch;
        } else
        {
            errPitches++;
        }

        summAMP +=val.amp;

        if (val.amp>maxAMP)
            maxAMP = val.amp;

        if (val.amp<minAMP)
            minAMP = val.amp;
    }

    float amp20perc = (maxAMP - minAMP)*0.2 + minAMP;
    float avgPitch = summPitch/(idx-errPitches);
    printf("AVG pitch %f amp %f [%f - %f] %f\n\n", avgPitch, summAMP/idx, minAMP, maxAMP, amp20perc );

    errPitches = 0;
    // find words and fill tempos
    std::vector<tempoidx> tempos;
    std::vector<float> wordPitches;
    std::vector<int> patternIndexes;

    uint patternLenMS = 0;
    uint fileRnd = 0;

/// inter {
    // pitch and amp interpolation
    for(idx=0; idx<pitchesSize-IFRAME_SIZE;idx++)
    {
        float spitch = 0, samp=0;
        for(int fdx=0; fdx<IFRAME_SIZE; fdx++)
        {
            amppitch val = pitches[idx+fdx];
            spitch += val.pitch;
            samp += val.amp;
        }
        pitches[idx].amp = samp/IFRAME_SIZE;
        pitches[idx].pitch = spitch/IFRAME_SIZE;
    }

/// inter }

    float wspitch=0;
    float wpitch = 0;
    beginOfWord = 0;
    for (idx=0; idx<peaks.size(); idx++)
    {
        endOfWord = deltaMS*peaks[idx]/1000;
        fileRnd = std::rand() % (patternFileNames.size());
        WavInFile patternfile(patternFileNames[fileRnd].c_str());
        patternLenMS = patternfile.getLengthMS();
        stretchFactor = patternLenMS/(endOfWord-beginOfWord)/1000.0f;
        tempoidx val;
        val.tempo = (stretchFactor-1.0)/0.01;

        tempos.push_back(val);
        patternIndexes.push_back(fileRnd);

        printf("[%d] Word at the [%.3f-%.3f] len: %.3f with %.3f\n", idx, beginOfWord, endOfWord, endOfWord-beginOfWord, stretchFactor);


        int sdx = 0;
        for (sdx=beginIdx; sdx<peaks[idx]; sdx++)
        {
            amppitch voicedatas = pitches[sdx];
            amppitch nextvoice = pitches[sdx+1];
            if (voicedatas.pitch != 0)
            {
                wspitch += voicedatas.pitch;
            } else
                errPitches++;
        }


        wpitch = wspitch/(sdx-beginIdx-errPitches);
        wordPitches.push_back(wpitch);

        wspitch = 0;
        errPitches = 0;
        //must add to vector
        beginOfWord = endOfWord;
        beginIdx = peaks[idx];
    }
    ///// end of searching words


    SAMPLETYPE *sampleBuffer;
    try 
    {
        int countSamples = 0, summSamples = 0;
        int sdx=0;
        int saveSamples = 0;
        bool rewind = false;

        std::cout << "tempos size " << tempos.size() << std::endl;
        for (ulong idx=0; idx<tempos.size(); idx++)
        {
            uint patternIndex = patternIndexes[idx];
            WavInFile patternfile(patternFileNames[patternIndex].c_str());
            uint numSamples = patternfile.getNumSamples();
            sampleBuffer = new float[numSamples];
            // Read a chunk of samples from the input file
            int nSamples = patternfile.read(sampleBuffer, numSamples);

            pitchFactor = 1;
            SoundTouch soundTouch;
            // Setup the 'SoundTouch' object for processing the sound
            soundTouch.setSampleRate((int)patternfile.getSampleRate());
            soundTouch.setChannels(1);

            tempoidx val = tempos[idx];
            float wpitch = wordPitches[idx];

            float factor = 1.0 + 0.01 * val.tempo;
            int outSize = ceil(numSamples/factor);

            std::cout << "here " << val.tempo << " " << wpitch << std::endl;
            // here processing
            SAMPLETYPE *outBuffer = new float[outSize];

            // we can change params
            soundTouch.setTempoChange(val.tempo);
            soundTouch.setPitch(1-(1-(wpitch/avgPitch))/prettifyPitch);

            // Feed the samples into SoundTouch processor
            soundTouch.putSamples(sampleBuffer, nSamples);


            delete[] sampleBuffer;

            // Now the input file is processed, yet 'flush' few last samples that are
            // hiding in the SoundTouch's internal processing pipeline.
            soundTouch.flush();

            int outSamples  = soundTouch.receiveSamples(outBuffer, outSize);
            /// pitching
            summSamples += outSamples;
            int ptr = 0;

            amppitch voicedatas = pitches[sdx];
            amppitch nextvoice;
            amppitch patterndatas;

            for (;;)
            {
                if (sdx<pitchesSize)
                    voicedatas = pitches[sdx];
                if (sdx+1<pitchesSize)
                    nextvoice = pitches[sdx+1];
                else
                    nextvoice = voicedatas;

                if (countSamples+bufferLengthFrames > summSamples)
                {
                    saveSamples = summSamples-countSamples;

                    rewind = true;
                    break; // out from here
                }

                if (rewind)
                {
                    saveSamples = bufferLengthFrames-saveSamples;
                    rewind = false;
                } else
                {
                    saveSamples = bufferLengthFrames;
                }
                countSamples+=bufferLengthFrames;

                // process data
                patterndatas = dywapitch_computepitch(&patterntracker, outBuffer+ptr, 0, saveSamples);
                // need to convert amp to factor. i.e. get main amp of cat's sound and find diff with original voice

                if (voicedatas.amp!=0 && patterndatas.amp!=0)
                {
                    ampFactor = 1-(1-(voicedatas.amp/patterndatas.amp))/prettifyAmp; // factor
                }

                if (ampFactor>0.9)
                    ampFactor = 0.9;

                // need to convert pitch to factor. i.e. get main pitch of cat's sound and find diff with original voice
                if (voicedatas.pitch!=0 && patterndatas.pitch!=0)
                {
                    pitchFactor = 1-(1-((voicedatas.pitch+prevdatas.pitch)/2)/patterndatas.pitch)/3; // factor

                    if (pitchFactor>1)
                        pitchFactor = 1;
                }

                float deltaAmp = (ampFactor-prevAmpFactor)/saveSamples;
                float curAmpFactor = prevAmpFactor;
                if (voicedatas.amp < 0.01)
                {
                    curAmpFactor = 0;
                    deltaAmp = 0;
                }

                for (int adx = 0; adx < saveSamples; adx++)
                {
                    outBuffer[ptr+adx] *= curAmpFactor;
                    curAmpFactor += deltaAmp;
                }

                prevAmpFactor = ampFactor;
                // end of processing


                // move pointer
                ptr+=saveSamples;
                sdx++;
                if (sdx == pitchesSize)
                    break;
            }

            outFile.write(outBuffer, outSamples);

            delete[] outBuffer;
            outBuffer = NULL;

        }
        printf("SAMPLES %d\n", sdx);
        delete[] voiceData;
    } 
    catch (const runtime_error &e) 
    {
        delete[] voiceData;
        // An exception occurred during processing, display an error message
        fprintf(stderr, "%s\n", e.what());
        return -1;
    }

    return 0;
}
