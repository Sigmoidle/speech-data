%Turn .wav to MFCC spectogram


names = ["albi","alejandro","alex",...
       "alexander","aurelie","benjamin",...
       "brennan", "felipe", "harry",...
       "hemal", "hugo", "max",...
       "nathaniel", "owen", "ruaridh",...
       "ruby", "sarah", "sophie",...
       "vav", "yan"];

melFilterBands = 40;

for i=1:20
    filename = strcat(names(i),".wav");
    name = names(i);
    [data, sampleRate] = audioread(filename, "double");
    
    totalSamples = audioinfo(filename).TotalSamples;
    samplesPerFrame = int16(totalSamples/((totalSamples/sampleRate)*10)); %exactly 10ms worth of samples
    hopSize = int16(samplesPerFrame/2);
    numFrames = floor((totalSamples/samplesPerFrame)*2)-2;
    
    processFrames(samplesPerFrame, hopSize, numFrames, data, melFilterBands, sampleRate, name)
end


function writeHTKFile(mfccSpectogram, numFrames, hopSize, sampleRate, filename)
    mfcFileName = strcat(filename, ".mfc");
    fid = fopen(mfcFileName,"w", "ieee-be");

    %write header
    fwrite(fid, numFrames, "int32");
    fwrite(fid, hopSize/sampleRate*1000000, "int32");
    fwrite(fid, 39*4, "int16");
    fwrite(fid, 6, "int16");

    for i=1: numFrames
        for j=1:39

            fwrite(fid, mfccSpectogram(j,i), "float32");
        end
    end
    fid = fclose(fid);
end

function mfccSpectogramWithTemporalInformation = processFrames(samplesPerFrame, hopSize, numFrames, data, bands, sampleRate, filename)
    mfccSpectogram = zeros(13, numFrames);
    for frameNumber = 1:numFrames
        dataSlice = data((frameNumber)*hopSize-hopSize+1:((frameNumber)*hopSize)+samplesPerFrame-hopSize+1);
        mfcc = makeMFCC(dataSlice, samplesPerFrame, bands, sampleRate);
        mfccSpectogram(:,frameNumber)=mfcc';
    end
    mfccSpectogramDelta = sgolayfilt(mfccSpectogram,1,9);
    mfccSpectogramDeltaDelta = sgolayfilt(mfccSpectogram,2,9);
    mfccSpectogramWithTemporalInformation = [mfccSpectogram;mfccSpectogramDelta;mfccSpectogramDeltaDelta];
    writeHTKFile(mfccSpectogramWithTemporalInformation, numFrames, hopSize, sampleRate, filename);
end


function mfcc = makeMFCC(data, samples, bands, sampleRate)
    %1st apply hamming window function
    dataHammed = hamming(samples)'.*data;
    %2nd apply fft to data and get mag phase 
    [mag, ~] = getMagPhase(dataHammed, samples);
    %3rd convert magnitude to decibels
    magDB = mag2db(mag);
    %4th run the magnitude vectors though the mel-filterbank
    filterBands = createMelFrequencyBands(sampleRate, bands);
    melFilterBank = applyMelFilter(magDB, filterBands, bands);
    %5th perform dct on mel filter bank
    dctOut = applyDCT(melFilterBank);
    %6th truncate
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

function bands = createMelFrequencyBands(sampleRate, noOfBands)
    highestFrequency = sampleRate/2;
    lowestMel = frequencyToMelScale(0);
    highestMel = frequencyToMelScale(highestFrequency);
    melBands = linspace(lowestMel, highestMel, noOfBands+2);
    frequencyBands = zeros(size(melBands));
    %convert mel bands back into frequency values
    for band = 1:noOfBands+2
        frequencyBands(band) = melScaleToFrequency(melBands(band));
    end
    bands = frequencyBands;
end

function sample = hzToSample(hz)
    sample = round(hz/10);
end

function melFilterBank = applyMelFilter(data, bands, noOfBands)
    melFilterBank = zeros([1,noOfBands]);
    for band = 2:noOfBands+1
        filter = zeros([1,800]);
        middlePoint = hzToSample(bands(band));
        backPoint = hzToSample(bands(band-1))+1;%+1 because matlab weird
        frontPoint = hzToSample(bands(band+1));
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
