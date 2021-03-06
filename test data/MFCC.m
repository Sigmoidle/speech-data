%Turn .wav to MFCC spectogram


names = ["albi","alejandro","alex",...
       "alexander","aurelie","benjamin",...
       "brennan", "felipe", "harry",...
       "hemal", "hugo", "max",...
       "nathaniel", "owen", "ruaridh",...
       "ruby", "sarah", "sophie",...
       "vav", "yan"];

%the number of mel bands to use, however the number of dct coefficients
%will still be 13 and there'll be 3 lots of them because we take two
%delters of first and second order to augiment the MFCC spectogram with
melFilterBands = 40;
frameSize = 20;

%This for loop iterates over every name in the names list, extracts
%important and specific information about the .wav file then calls for them
%to be processed and saved into MFCC spectogram HTK files
for i=1:20
    filename = strcat(names(i),".wav");
    name = names(i);
    [data, sampleRate] = audioread(filename, "double");
    
    totalSamples = audioinfo(filename).TotalSamples;
    samplesPerFrame = int16(totalSamples/((totalSamples/sampleRate)*5)); %exactly 20ms worth of samples for 16000hz sample rate
    hopSize = 10;
    numFrames = floor((totalSamples/sampleRate*1000)/frameSize)-2;
    
    processFrames(samplesPerFrame, hopSize, numFrames, data, melFilterBands, sampleRate, name)
end

%This function writes the NAME.mfc HTF format file
function writeHTKFile(mfccSpectogram, numFrames, hopSize, sampleRate, filename)
    mfcFileName = strcat(filename, ".mfc");
    fid = fopen(mfcFileName,"w", "ieee-be");

    %write header
    fwrite(fid, int32(numFrames), "int32");
    fwrite(fid, 100000, "int32");
    fwrite(fid, 39*4, "int16");
    fwrite(fid, int16(6+256+512), "int16");
    
    %write individual lines of data
    for i=1: numFrames
        for j=1:39
            fwrite(fid, mfccSpectogram(j,i), "float32");
        end
    end
    fid = fclose(fid);
end

%"ProcessFrames"
%This function goes through the .wav file and calls the "makeMFCC" function
%on each "Frame" the frames are overlapping too. Each time the MFCC is
%calculated for a fame it will add it to the overall MFCC spectogram that
%is being made.
%After the full MFCC spectogram is made for the .wav file, it's first
%and second temporal derivitives are augimented onto it to add additional
%information for the HMM
%Lastly "writeHTKFile" is called to write the file
function mfccSpectogramWithTemporalInformation = processFrames(samplesPerFrame, hopSize, numFrames, data, bands, sampleRate, filename)
    mfccSpectogram = zeros(13, numFrames);
    for frameNumber = 1:numFrames
        dataSlice = data((frameNumber)*hopSize-hopSize+1:((frameNumber)*hopSize)+samplesPerFrame-hopSize+1);
        mfcc = makeMFCC(dataSlice, samplesPerFrame, bands, sampleRate);
        mfccSpectogram(:,frameNumber)=mfcc';
    end
    %Apply a sgolayfilter to the MFCC spectogram to get delta and delta
    %delta information from the data. This is then augimented onto the
    %original MFCC spectogram
    mfccSpectogramDelta = sgolayfilt(mfccSpectogram,1,9);
    mfccSpectogramDeltaDelta = sgolayfilt(mfccSpectogram,2,9);
    mfccSpectogramWithTemporalInformation = [mfccSpectogram;mfccSpectogramDelta;mfccSpectogramDeltaDelta];
    mfccSpectogramWithTemporalInformation(isnan(mfccSpectogramWithTemporalInformation)) = 0;
    writeHTKFile(mfccSpectogramWithTemporalInformation, numFrames, hopSize, sampleRate, filename);
end

%Turns a frame of time domain .wav file into Xms of MFCC data
function mfcc = makeMFCC(data, samples, bands, sampleRate)
    %1st apply hamming window function
    dataHammed = hamming(samples)'.*data;
    %2nd apply fft to data and get mag phase (discard phase as it's not
    %needed)
    [mag, ~] = getMagPhase(dataHammed, samples);
    %3rd convert magnitude to decibels
    magDB = mag2db(mag);
    %4th run the magnitude vectors though the mel-filterbank
    filterBands = createMelFrequencyBands(sampleRate, bands);
    melFilterBank = applyMelFilter(magDB, filterBands, bands, samples);
    %5th perform dct on mel filter bank, this is to make useful information
    %about tembra etc very pronounced and less useful information like
    %volume not appear in the data.
    dctOut = applyDCT(melFilterBank);
    %6th truncate (remove the first co-efficient and take the next 13
    %co-efficiencts)
    mfcc = truncate(dctOut);
end

function [mag, phase] = getMagPhase(data, samplesPerFrame)
    dataDFT = fft(data);
    mag = abs(dataDFT(1:samplesPerFrame/2));
    phase = angle(dataDFT); 
end

function mel = frequencyToMelScale(f)
    mel = 2595*log10(1+(f/700));
end

function f = melScaleToFrequency(mel)
    f = 700*(power(10,mel/2595)-1);
end

%This function creates the mel-scaled frequency bands required for making a
%mel-scaled frequency bank 
function bands = createMelFrequencyBands(sampleRate, noOfBands)
    highestFrequency = sampleRate/2;
    lowestMel = frequencyToMelScale(0);
    highestMel = frequencyToMelScale(highestFrequency);
    %Using the highest and lowest mel-scailed frequencies the next line of
    %code makes a linearly interpolated bunch of points between those two
    %values (plus two for the start and the end)
    melBands = linspace(lowestMel, highestMel, noOfBands+2);
    frequencyBands = zeros(size(melBands));
    %converts mel bands back into frequency values
    for band = 1:noOfBands+2
        frequencyBands(band) = melScaleToFrequency(melBands(band));
    end
    bands = frequencyBands;
end

function sample = hzToSample(hz, samples)
    sample = round(hz/(8000/(samples/2)));
end

%This function uses the mel-scale frequency bands to create multiple
%mel-scaled filter banks (triangular overlapping). 
function melFilterBank = applyMelFilter(data, bands, noOfBands, samples)
    melFilterBank = zeros([1,noOfBands]);
    for band = 2:noOfBands+1
        filter = zeros([1,samples/2]);
        middlePoint = hzToSample(bands(band), samples);
        backPoint = hzToSample(bands(band-1), samples)+1;%+1 because matlab starts indexes at 1
        frontPoint = hzToSample(bands(band+1), samples);
        filter(backPoint:middlePoint) = linspace(0, 1, middlePoint-backPoint+1);
        filter(middlePoint:frontPoint) = linspace(1, 0, frontPoint-middlePoint+1);
        melFilterBank(band-1) = sum(data.*filter, "all");
    end
end

function dctOut = applyDCT(data)
    dctOut = dct(data);
end

function actualMFCC = truncate(data)
    actualMFCC = data(2:14);
end
