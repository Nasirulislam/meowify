/////////////////////////////////////////////////////////////////////
//
// The program reads a .wav sound
// file with mono-16bit-44100Hz sample format, process it
// and writes output into another .wav file.
//
/////////////////////////////////////////////////////////////////////

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include <cfloat>

#include "dywapitchtrack.h"
#include <stdexcept>
#include <vector>
#include <iostream>
#include "wavfile.h"

using namespace std;

// Time scaling factor, values > 1.0 increase, values < 1.0 decrease tempo
#define TIME_SCALE      0.85   // 15% slower tempo
// Processing sequence size (100 msec with 44100Hz samplerate)
#define SEQUENCE        4410
// Overlapping size (20 msec)
#define OVERLAP         882
// Best overlap offset seeking window (15 msec)
#define SEEK_WINDOW     662
// Processing sequence flat mid-section duration
#define FLAT_DURATION   (SEQUENCE - 2 * (OVERLAP))
// Theoretical interval between the processing seqeuences
#define SEQUENCE_SKIP   ((int)((SEQUENCE - OVERLAP) * (TIME_SCALE)))

typedef short SAMPLE;   // sample type, 16bit signed integer

// Use cross-correlation function to find best overlapping offset
// where input_prev and input_new match best with each other
int seek_best_overlap(const SAMPLE *input_prev, const SAMPLE *input_new)
{
   int i;
   int bestoffset = 0;
   float bestcorr = -1e30f;
   float temp[OVERLAP];

   // Precalculate overlapping slopes with input_prev
   for (i = 0; i < OVERLAP; i ++)
   {
      temp[i] = (float)(input_prev[i] * i * (OVERLAP - i));
   }

   // Find best overlap offset within [0..SEEK_WINDOW]
   for (i = 0; i < SEEK_WINDOW; i ++)
   {
      int j;
      float crosscorr = 0;

      for (j = 0; j < OVERLAP; j ++)
      {
         crosscorr += (float)input_new[i + j] * temp[j];
      }
      if (crosscorr > bestcorr)
      {
         // found new best offset candidate
         bestcorr = crosscorr;
         bestoffset = i;
      }
   }
   return bestoffset;
}


// Overlap 'input_prev' with 'input_new' by sliding the amplitudes during
// OVERLAP samples. Store result to 'output'.
void overlap(SAMPLE *output, const SAMPLE *input_prev, const SAMPLE *input_new)
{
   int i;

   for (i = 0; i < OVERLAP; i ++)
   {
      output[i] = (input_prev[i] * (OVERLAP - i) + input_new[i] * i) / OVERLAP;
   }
}


// SOLA algorithm. Performs time scaling for sample data given in 'input',
// write result to 'output'. Return number of output samples.
int sola(SAMPLE *output, const SAMPLE *input, int num_in_samples)
{
   int num_out_samples = 0;
   const SAMPLE *seq_offset = input;
   const SAMPLE *prev_offset;

   while (num_in_samples > SEQUENCE_SKIP + SEEK_WINDOW)
   {
      // copy flat mid-sequence from current processing sequence to output
      memcpy(output, seq_offset, FLAT_DURATION * sizeof(SAMPLE));
      // calculate a pointer to overlap at end of the processing sequence
      prev_offset = seq_offset + FLAT_DURATION;

      // update input pointer to theoretical next processing sequence begin
      input += SEQUENCE_SKIP - OVERLAP;
      // seek actual best matching offset using cross-correlation
      seq_offset = input + seek_best_overlap(prev_offset, input);

      // do overlapping between previous & new sequence, copy result to output
      overlap(output + FLAT_DURATION, prev_offset, seq_offset);

      // Update input & sequence pointers by overlapping amount
      seq_offset += OVERLAP;
      input  += OVERLAP;

      // Update output pointer & sample counters
      output += SEQUENCE - OVERLAP;
      num_out_samples += SEQUENCE - OVERLAP;
      num_in_samples -= SEQUENCE_SKIP;
   }

   return num_out_samples;
}



// Buffers for input/output sample data. For sake of simplicity, these are
// just made 'big enough' for the example purpose.
SAMPLE inbuffer[10240000];
SAMPLE outbuffer[20240000];

int main2(int numstr, char **pstr)
{

   if (numstr < 3)
   {
      printf("usage: solatest input.wav output.wav\n");
      return -1;
   }

   try
   {
      int insamples, outsamples;

      // Open input file
      WavInFile infile(pstr[1]);

      if ((infile.getSampleRate() != 44100) || (infile.getNumChannels() != 1))
      {
         printf("Sorry, this example processes mono audio sampled at 44100Hz.\n");
         return -1;
      }

      // Read data from input file
      insamples = infile.read(inbuffer, 10240000);

      // Process
      outsamples = sola(outbuffer, inbuffer, insamples);
      printf("Data to transform original %d %d\n", insamples, outsamples);
      // Write result to output file
      WavOutFile outfile(pstr[2], infile.getSampleRate(), infile.getNumBits(), infile.getNumChannels());
      outfile.write(outbuffer, outsamples);
   }
   catch (exception &e)
   {
      printf("Error: %s\n", e.what());
   }

   return 0;
}


void PitchShift(float pitchShift, long numSampsToProcess, long fftFrameSize, long osamp, float sampleRate, float *indata, float *outdata);

int main2(void)
{

    long numChannels = 1;
    long bufferLengthFrames = 512;
    dywapitchtracker voicetracker;
    dywapitch_inittracking(&voicetracker);

    dywapitchtracker patterntracker;
    dywapitch_inittracking(&patterntracker);

    float *voiceData = new float[bufferLengthFrames];
    float *patternData = new float[bufferLengthFrames];

    float semitones = 3.8;							// shift up by 3 semitones
    float pitchShift = pow(2.1f, semitones/12.0f);	// convert semitones to factor
    pitchShift = 1.377261;

    char voiceFileName[] = "input.wav";
    char patternFileName[] = "pattern.wav";
    char outFileName[] = "output.wav";

    // Create output file
    printf("Data to transform original sound\n");

    int idx = 0;
    float pitchFactor = 0;
    float ampFactor = 0;
    float prevAmpFactor = 0;
    int insamples, patternsamples;

    // Write result to output file
    WavOutFile outfile(outFileName, 44100, 16, numChannels);

    // Open input file
    WavInFile infile(voiceFileName);

    if ((infile.getSampleRate() != 44100) || (infile.getNumChannels() != 1))
    {
        printf("Sorry, this example processes mono audio sampled at 44100Hz for voice.wav.\n");
        return -1;
    }

    // Open parrern file
    WavInFile patternfile(patternFileName);

    if ((patternfile.getSampleRate() != 44100) || (patternfile.getNumChannels() != 1))
    {
        printf("Sorry, this example processes mono audio sampled at 44100Hz for meow.wav.\n");
        return -1;
    }

    std::vector<amppitch> pitches;

    for(;;)
    {
        // Read data from input file
        insamples = infile.read(voiceData, bufferLengthFrames);
        if (infile.eof())
        {
            break;
        }
        amppitch voicedatas = dywapitch_computepitch(&voicetracker, voiceData, 0, bufferLengthFrames);
        pitches.push_back(voicedatas);
    }
    infile.rewind();

    uint lenMS = infile.getLengthMS();
    float deltaMS = lenMS/float(pitches.size());
    std::cout << "pitches stores " << int(pitches.size()) << " numbers. " << lenMS << " " << deltaMS <<"\n";
    amppitch prevdatas;

    insamples = 0;

    float minAMP = FLT_MAX;
    float maxAMP = 0;
    float summAMP = 0;
    float summPitch = 0;
    float summPatternPitch = 0;

    int errPitches = 0;
    int errPatternPitches = 0;
    for(idx=0; idx<int(pitches.size());idx++)
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
    printf("AVG pitch %f amp %f [%f - %f] %f\n", avgPitch, summAMP/idx, minAMP, maxAMP, amp20perc );


    std::vector<amppitch> patternPitches;

    for(;;)
    {
        // Read data from pattern file
        patternsamples = patternfile.read(patternData, bufferLengthFrames);
        if (patternfile.eof())
        {

            break;
        }
        amppitch patterndatas = dywapitch_computepitch(&patterntracker, patternData, 0, bufferLengthFrames);
        patternPitches.push_back(patterndatas);
    }
    patternfile.rewind();

    uint patternLenMS = patternfile.getLengthMS();
    float patternDeltaMS = patternLenMS/float(patternPitches.size());
    int patternSize = int(patternPitches.size());
    std::cout << "pattern pitches stores " << patternSize << " numbers. " << patternLenMS << " " << patternDeltaMS <<"\n";

    for(idx=0; idx<int(patternPitches.size());idx++)
    {
        amppitch val = patternPitches[idx];
        if (val.pitch != 0)
        {
            summPatternPitch +=val.pitch;
        } else
        {
            errPatternPitches++;
        }
    }


    printf("AVG pattern pitch %f\n", /*avgPitch/*/(summPatternPitch/(idx-errPatternPitches)));

    float beginOfWord = 0;
    float endOfWord = 0;
    int beginIdx = 0;
    float stretchFactor = 0;
    for(idx=0; idx<int(pitches.size());idx++)
    {
        amppitch voicedatas = pitches[idx];
        amppitch nextvoice = pitches[idx+1];

        if (prevdatas.pitch == 0 && voicedatas.pitch != 0 && nextvoice.pitch != 0 && voicedatas.amp >= amp20perc)
        {
            endOfWord = deltaMS*idx/1000;
            stretchFactor = patternLenMS/(endOfWord-beginOfWord)/1000.0f;
            printf("[%d] Word at the [%.3f-%.3f] len: %.3f [%d] with %.3f \n", idx, beginOfWord, endOfWord, endOfWord-beginOfWord, idx-beginIdx, stretchFactor);
            //must add to vector
            beginOfWord = endOfWord;
            beginIdx = idx;
        }
        prevdatas = voicedatas;
    }
    endOfWord = lenMS/1000.0f;
    stretchFactor = patternLenMS/(endOfWord-beginOfWord)/1000.0f;
    printf("[%d] Word at the [%.3f-%.3f] len: %.3f [%d] with %.3f \n", idx, beginOfWord, endOfWord, endOfWord-beginOfWord, idx-beginIdx, stretchFactor);

    for(idx=0; idx<int(pitches.size());idx++)
    {
        amppitch voicedatas = pitches[idx];
        amppitch nextvoice = pitches[idx+1];
        if (prevdatas.pitch == 0 && voicedatas.pitch != 0 && nextvoice.pitch != 0 && voicedatas.amp >= amp20perc)
        {
            patternfile.rewind();
        }

        // Read data from pattern file
        patternsamples = patternfile.read(patternData, bufferLengthFrames);
        if (patternfile.eof())
        {
            patternfile.rewind();
            infile.read(patternData+patternsamples, bufferLengthFrames-patternsamples);
        }

        // computes the pitch. Pass the inited dywapitchtracker structure
        // samples : a pointer to the sample buffer
        // startsample : the index of teh first sample to use in teh sample buffer
        // samplecount : the number of samples to use to compte the pitch
        // return 0.0 if no pitch was found (sound too low, noise, etc..)
        // double dywapitch_computepitch(dywapitchtracker *pitchtracker, double * samples, int startsample, int samplecount);

        amppitch patterndatas = dywapitch_computepitch(&patterntracker, patternData, 0, bufferLengthFrames);

        // need to convert amp to factor. i.e. get main amp of cat's sound and find diff with original voice
        if (voicedatas.amp!=0 && patterndatas.amp!=0)
        {
            ampFactor = voicedatas.amp/patterndatas.amp; // factor
        }

        if (ampFactor>0.9)
            ampFactor = 0.9;

        // need to convert pitch to factor. i.e. get main pitch of cat's sound and find diff with original voice
        if (voicedatas.pitch!=0 && patterndatas.pitch!=0)
        {
            pitchFactor = (voicedatas.pitch+prevdatas.pitch)/2 /patterndatas.pitch; // factor
//            pitchFactor = voicedatas.pitch/avgPitch;

            if (pitchFactor>1)
                pitchFactor = 1;
        }

        float delta = (1 - pitchFactor);
//        printf("factor is %f -> %f = %f \n", prevAmpFactor, ampFactor, (ampFactor-prevAmpFactor)/bufferLengthFrames);

//        printf("pitch is %f, pitch is %f factor must be %f ampFactor=%f\n", voicedatas.pitch, patterndatas.pitch, pitchFactor, ampFactor);
//        printf("factor must be %f ampfactor=%f delta=%f (%f)\n", pitchFactor, ampFactor, delta, 1-delta*ampFactor);

        // void PitchShift(float pitchShift, long numSampsToProcess, long fftFrameSize, long osamp, float sampleRate, float *indata, float *outdata)
        // --------------------------------- call PitchShift() ---------------------------------
//        PitchShift(1-(1-pitchFactor)/3, bufferLengthFrames, bufferLengthFrames, 4, mAiffGetSampleRate(patternFileName), patternData, patternData);
//        PitchShift(pitchFactor, bufferLengthFrames, bufferLengthFrames, 4, mAiffGetSampleRate(patternFileName), patternData, patternData);

        PitchShift(1-(delta)/2.5, bufferLengthFrames, bufferLengthFrames, 4, 44100, patternData, patternData);

        // ----------------------------------------------------------------------------------------
        float deltaAmp = (ampFactor-prevAmpFactor)/bufferLengthFrames;
        float curAmpFactor = prevAmpFactor;

        for (int idx = 0; idx < bufferLengthFrames; idx++)
        {
            patternData[idx] *= curAmpFactor;
            curAmpFactor += deltaAmp;
        }

        prevAmpFactor = ampFactor;
        prevdatas = voicedatas;

        // write processed data to output file
        outfile.write(patternData, bufferLengthFrames);
    }

    delete[] voiceData;
    delete[] patternData;

    // That's all there is to it.
    printf("You have %d frames\n", idx);

//    system("open out.aif");


    return 0;
}
