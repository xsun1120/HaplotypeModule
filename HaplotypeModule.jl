#Pkg.clone("https://github.com/JuliaLang/DataStructures.jl.git")

module HaplotypeModule

#using DataStructures.OrderedDict

type SNPInfo
    name::ASCIIString
    chr::Int64
    pos::Int64
    index::Int64
    seg::Int64
    segOnChr::Int64
    freq::Float64
    breedFreq::Array{Float64}
end

type HapInfo
    id::ASCIIString
    index::Int64
    alleleState::ASCIIString
    #seg::Int64
    #segOnChr::Int64
    count::Int64
    freq::Float64
    rare::Bool
    breedFreq::Array{Float64}
end

type SegInfo
    index::Int64
    chr::Int64
    indexOnChr::Int64
    numberSNPs::Int64
    numberUnique::Int64
    numberCommon::Int64
    begSNP::ASCIIString
    endSNP::ASCIIString
    length::Int64
    uniqueHap::Dict{ASCIIString, HapInfo}
    commonHap::Set{ASCIIString}
end

type Individual
    name::ASCIIString
    breed::Int64
    index::Int64
    matHap::ASCIIString
    patHap::ASCIIString
    geno::Array{Float64}
    haplo::Array{Float64}
    compSNP::Array{Float64}
    compHap::Array{Float64}
end

type Genome
    numberBreeds::Int64
    Breeds::Array{ASCIIString}
    BreedSizes::Array{Int64}
    numberSNPs::Int64
    Chromosomes::Array{Int64}
    SNPInfoMap::Dict{ASCIIString, SNPInfo}
    SNPInfoVec::Array{ASCIIString}
    Segments::Array{SegInfo}
    totalNumberUnique::Int64
    totalNumberCommon::Int64
    Pop::Dict{ASCIIString, Individual}
    numberInd::Int64
end

function makeGenome(breedStr::ASCIIString)  # make a genome object
    # breedStr: 3-letter abbr. of breed names to include in analysis, separated by , e.g. "AAN,CHA,GVH"
    breedData = split(breedStr, ',')
    nBreeds = length(breedData)
    breeds = ASCIIString[]
    for br = 1:nBreeds
        push!(breeds, ascii(breedData[br]))
    end
    println("$nBreeds breeds: $breeds")

    numberSNPs = 0::Int64
    Chromosomes = Int64[]
    SNPInfoMap = Dict{ASCIIString, SNPInfo}()
    SNPInfoVec = ASCIIString[]
    Segments = SegInfo[]
    totalNumberUnique = 0::Int64
    totalNumberCommon = 0::Int64
    Pop = Dict{ASCIIString, Individual}()
    numberInd = 0::Int64

    genome = Genome(nBreeds, breeds, zeros(Int64,nBreeds), numberSNPs, Chromosomes, SNPInfoMap, SNPInfoVec, Segments,
    totalNumberUnique, totalNumberCommon, Pop, numberInd)
    
    return genome
end


function inputPhenotype!(g::Genome, fileName::ASCIIString, breedIdx::Int64)
    # input ids for animals in the breedIdx-th breed from fileName in GenSel phenotype file format
    d = open(fileName)
    lineNumber = 0
    for ln in eachline(d)
        lineNumber += 1
        (lineNumber == 1) && continue
        lnData = split(ln)
        tempInd = ascii(lnData[1])
        if !haskey(g.Pop, tempInd)
            g.Pop[tempInd] = Individual(tempInd, breedIdx, 0, "", "", [], [], [], [])
        else
            error("Individual $tempInd already exists!")
        end
    end
    close(d)
    g.BreedSizes[breedIdx] = lineNumber-1
    g.numberInd += g.BreedSizes[breedIdx]
    println("Number of individuals in breed $(g.Breeds[breedIdx]) : $(g.BreedSizes[breedIdx])")
end

function inputSNPInfo!(g::Genome, fileName::ASCIIString)
    # input sorted SNP id, chr and pos information from map file fileName
    # these SNPs are used to construct haplotypes and their genotypes should be available in all breeds
    d = open(fileName)
    lineNumber = 0
    thisChr = 0
    snpIdx = 0
    chrIdx = 0
    posIdx = 0
    for ln in eachline(d)
        lineNumber += 1
        lnData = split(ln)
        if lineNumber == 1
            snpIdx = find(x->x=="SNP",lnData)[1]
            chrIdx = find(x->x=="Chromosome",lnData)[1]
            posIdx = find(x->x=="Position_bp",lnData)[1]
            continue
        end

        tempSNP = ascii(lnData[snpIdx])
        tempChr = parse(Int64, lnData[chrIdx])
        tempPos = parse(Int64, lnData[posIdx])
        if tempChr != thisChr
            push!(g.Chromosomes, tempChr)
            thisChr = tempChr
        end
        
        if !haskey(g.SNPInfoMap,tempSNP)
            g.SNPInfoMap[tempSNP] = SNPInfo(tempSNP, tempChr, tempPos, lineNumber-1, 0, 0, 0.0, zeros(g.numberBreeds))
            push!(g.SNPInfoVec, tempSNP)
        else
            error("SNP $tempSNP duplicated!")
        end
    end
    close(d)
    g.numberSNPs = length(g.SNPInfoMap)
    println("Number of chromosomes: $(length(g.Chromosomes))")
    println("Read $(g.numberSNPs) SNPs from $fileName")
end

function genomeSegmentation!(g::Genome, segLength::Int64, fileName::ASCIIString)
    # divide SNPs into windows with segLength basepairs, results written to fileName.* 
    numberSegments = 0;
    numberSegmentsThisChr = 0;
    println("Genome divided into $segLength bp segments.")
    outSNP = open("$fileName.seg","w")
    write(outSNP, "SNP chr pos index seg segOnChr\n")
    tempChr = 0
    thisSeg = 0
    prevSeg = 0
    for tempSNP in g.SNPInfoVec
        if g.SNPInfoMap[tempSNP].chr != tempChr
            tempChr = g.SNPInfoMap[tempSNP].chr
            #println("Number of segments on chromosome $tempChr: $numberSegmentsThisChr")
            numberSegments += numberSegmentsThisChr
            numberSegmentsThisChr = 0
            thisSeg = 0
            prevSeg = 0
        end

        tempSeg = ceil(Int64, g.SNPInfoMap[tempSNP].pos/segLength)
        if tempSeg != prevSeg
            prevSeg = tempSeg
            thisSeg += 1
            numberSegmentsThisChr += 1
        end
        g.SNPInfoMap[tempSNP].segOnChr = thisSeg
        g.SNPInfoMap[tempSNP].seg = thisSeg + numberSegments
        write(outSNP, join((g.SNPInfoMap[tempSNP].name, g.SNPInfoMap[tempSNP].chr, g.SNPInfoMap[tempSNP].pos,
              g.SNPInfoMap[tempSNP].index, g.SNPInfoMap[tempSNP].seg, g.SNPInfoMap[tempSNP].segOnChr), " "),"\n")
    end
    numberSegments += numberSegmentsThisChr
    close(outSNP)

    # fill Segments
    outSeg = open("$fileName.ss","w")
    thisSeg = 0
    numberSNPsThisSeg = 0
    for j = 1:length(g.SNPInfoVec)
        tempSeg = g.SNPInfoMap[g.SNPInfoVec[j]].seg
        if tempSeg != thisSeg
            tempSegInfo = SegInfo(tempSeg, g.SNPInfoMap[g.SNPInfoVec[j]].chr, g.SNPInfoMap[g.SNPInfoVec[j]].segOnChr,
            0, 0, 0, "", "", 0, Dict{ASCIIString, HapInfo}(), Set{ASCIIString}())
            push!(g.Segments, tempSegInfo)
            g.Segments[end].begSNP = g.SNPInfoVec[j]
            if thisSeg > 0
                g.Segments[end-1].endSNP = g.SNPInfoVec[j-1]
                g.Segments[end-1].length = g.SNPInfoMap[g.Segments[end-1].endSNP].pos - g.SNPInfoMap[g.Segments[end-1].begSNP].pos + 1
                g.Segments[end-1].numberSNPs = numberSNPsThisSeg
            end
            thisSeg = tempSeg
            numberSNPsThisSeg = 1
        else
            numberSNPsThisSeg += 1
        end
    end
    g.Segments[end].endSNP = g.SNPInfoVec[end]
    g.Segments[end].length = g.SNPInfoMap[g.Segments[end].endSNP].pos - g.SNPInfoMap[g.Segments[end].begSNP].pos + 1
    g.Segments[end].numberSNPs = numberSNPsThisSeg
    println("Size of Segments: ",length(g.Segments))

    println("Total number of segments: $numberSegments")
    write(outSeg, "index indexOnChr chr nSNPs begPos endPos length\n")
    for tempSegInfo in g.Segments
        write(outSeg, join((tempSegInfo.index, tempSegInfo.indexOnChr, tempSegInfo.chr, tempSegInfo.numberSNPs,
        g.SNPInfoMap[tempSegInfo.begSNP].pos, g.SNPInfoMap[tempSegInfo.endSNP].pos, tempSegInfo.length), " "),"\n")
    end
    close(outSeg)
end


function inputPhasedGenotype!(g::Genome, breedIdx::Int64, fileName::ASCIIString)
    # input phased genotypes for individuals in the breedIdx-th breed from fileName
    # this function reads *.vcf.gz format for each chromosome
    indColumnIdx = []
    indNames = []
    run(`gzip -d $fileName.gz`)
    d = open(fileName)
    lineNumber = 0
    for ln in eachline(d)
        lineNumber += 1
        (lineNumber < 10) && continue
        if lineNumber == 10
            lnData = split(ln)
            for i = 10:length(lnData)
                tempInd = ascii(lnData[i])
                haskey(g.Pop, tempInd) || continue
                # g.Pop[tempInd].index = i
                push!(indColumnIdx, i)
                push!(indNames, tempInd)
            end
            println("Number of individuals found in $fileName: $(length(indColumnIdx))")
            continue
        end

        lnData = split(ln)
        tempSNP = ascii(lnData[3])
        haskey(g.SNPInfoMap, tempSNP) || continue
        #for tempInd in keys(g.Pop)
            #g.Pop[tempInd].breed == breedIdx || continue
            #idx = g.Pop[tempInd].index
        for i = 1:length(indColumnIdx)
            tempInd = indNames[i]
            idx = indColumnIdx[i]

            a1 = ascii(lnData[idx])[1]
            a2 = ascii(lnData[idx])[3]
            g.Pop[tempInd].matHap *= string(a1) # join([g.Pop[tempInd].matHap, a1])
            g.Pop[tempInd].patHap *= string(a2) # join([g.Pop[tempInd].patHap, a2])

            g.SNPInfoMap[tempSNP].breedFreq[breedIdx] += (parse(Int64, a1) + parse(Int64, a2))
            g.SNPInfoMap[tempSNP].freq += (parse(Int64, a1) + parse(Int64, a2))
        end
    end
    run(`gzip $fileName`)
    println("Read phased genotypes for $(lineNumber-10) SNPs from $fileName")
    close(d)    
end

function getHaplotypeAlleles!(g::Genome, rareHapFreq::Float64, fileName::ASCIIString)
    # construct haplotype alleles with length previously specified
    # common haplotypes have frequency larger than rareHapFreq
    # results written to fileName.*
    outSeg = open("$fileName.unique","w")
    write(outSeg, "index chr nUnique nCommon\n")
    for j = 1:length(g.Segments)
        begIdx = g.SNPInfoMap[g.Segments[j].begSNP].index
        endIdx = g.SNPInfoMap[g.Segments[j].endSNP].index
        for tempInd in keys(g.Pop)
            br = g.Pop[tempInd].breed
            # maternal
            tempHap = g.Pop[tempInd].matHap[begIdx:endIdx]
            if !haskey(g.Segments[j].uniqueHap, tempHap)
                tempHapInfo = HapInfo("", 0, tempHap, 1, 0.0, false, zeros(g.numberBreeds))
                tempHapInfo.breedFreq[br] = 1
                g.Segments[j].uniqueHap[tempHap] = tempHapInfo
            else
                g.Segments[j].uniqueHap[tempHap].count += 1
                g.Segments[j].uniqueHap[tempHap].breedFreq[br] += 1
            end
            # paternal
            tempHap = g.Pop[tempInd].patHap[begIdx:endIdx]
            if !haskey(g.Segments[j].uniqueHap, tempHap)
                tempHapInfo = HapInfo("", 0, tempHap, 1, 0.0, false, zeros(g.numberBreeds))
                tempHapInfo.breedFreq[br] = 1
                g.Segments[j].uniqueHap[tempHap] = tempHapInfo
            else
                g.Segments[j].uniqueHap[tempHap].count += 1
                g.Segments[j].uniqueHap[tempHap].breedFreq[br] += 1
            end
        end

        g.Segments[j].numberUnique = length(g.Segments[j].uniqueHap)
        for tempHap in keys(g.Segments[j].uniqueHap)
            g.Segments[j].uniqueHap[tempHap].freq = g.Segments[j].uniqueHap[tempHap].count / (2.0*g.numberInd)
            g.Segments[j].uniqueHap[tempHap].breedFreq ./= (2.0*g.BreedSizes)
            g.Segments[j].uniqueHap[tempHap].rare = g.Segments[j].uniqueHap[tempHap].freq <= rareHapFreq
            g.Segments[j].uniqueHap[tempHap].rare || push!(g.Segments[j].commonHap, tempHap)
        end
        g.Segments[j].numberCommon = length(g.Segments[j].commonHap)

        write(outSeg, join((g.Segments[j].index, g.Segments[j].chr, 
        g.Segments[j].numberUnique, g.Segments[j].numberCommon), " "),"\n")

        g.totalNumberUnique += g.Segments[j].numberUnique
        g.totalNumberCommon += g.Segments[j].numberCommon
    end
    close(outSeg)
    println("Total number of unique haplotypes: $(g.totalNumberUnique)")
    println("Total number of common haplotypes: $(g.totalNumberCommon)")
end


function getBreedComposition!(g::Genome, hapLength::Int64, rareHapFreq::Float64, mapFile::ASCIIString, outFile::ASCIIString)
    # calculate breed composition using both SNPs and haplotypes, calling all above functions
    inputSNPInfo!(g, mapFile)
    genomeSegmentation!(g, hapLength, outFile)

    for br = 1:(g.numberBreeds)
        tempBreed = g.Breeds[br]
        phenFile = "$tempBreed.all50k"
        println("Read individuals of breed $tempBreed from $phenFile")
        inputPhenotype!(g, phenFile, br)
        for chr = 1:length(g.Chromosomes)
            genoFile = "phased.HD.ab-genotype-HD-$tempBreed.$(g.Chromosomes[chr]).vcf.small"
            inputPhasedGenotype!(g, br, genoFile)
        end
    end
    println("Total number of individuals: $(g.numberInd)")
    println("Size of Pop: $(length(g.Pop))")

    getHaplotypeAlleles!(g, rareHapFreq, outFile)

    outHap = open("$outFile.freqHap","w")
    write(outHap, "index seg chr segChr freq $(join(g.Breeds," ")) states\n")

    XHap = zeros(Float64, g.totalNumberCommon, g.numberBreeds)
    rowIdx = 0
    for j = 1:length(g.Segments)
        for tempHap in g.Segments[j].commonHap
            rowIdx += 1
            g.Segments[j].uniqueHap[tempHap].index = rowIdx
            XHap[rowIdx,:] = g.Segments[j].uniqueHap[tempHap].breedFreq
            write(outHap, "$rowIdx $(g.Segments[j].index) $(g.Segments[j].chr) $(g.Segments[j].indexOnChr) $(g.Segments[j].uniqueHap[tempHap].freq) $(join(g.Segments[j].uniqueHap[tempHap].breedFreq," ")) $(g.Segments[j].uniqueHap[tempHap].alleleState)\n" )
        end
    end
    close(outHap)
    println("Haplotype frequencies written to $outFile.freqHap")

    outSNP = open("$outFile.freqSNP","w")
    write(outSNP, "SNP chr pos freq $(join(g.Breeds," "))\n")
    XSNP = zeros(Float64, g.numberSNPs, g.numberBreeds)
    for tempSNP in g.SNPInfoVec
        g.SNPInfoMap[tempSNP].freq /= (2.0*g.numberInd)
        g.SNPInfoMap[tempSNP].breedFreq ./= (2.0*g.BreedSizes)
        write(outSNP, "$tempSNP $(g.SNPInfoMap[tempSNP].chr) $(g.SNPInfoMap[tempSNP].pos) $(g.SNPInfoMap[tempSNP].freq) $(join(g.SNPInfoMap[tempSNP].breedFreq," "))\n")
        rowIdx = g.SNPInfoMap[tempSNP].index
        XSNP[rowIdx,:] = g.SNPInfoMap[tempSNP].breedFreq
    end
    close(outSNP)
    println("SNP frequencies written to $outFile.freqSNP")
    
    outHap = open("$outFile.compHap","w")
    outSNP = open("$outFile.compSNP","w")
    write(outHap, "id breed $(join(g.Breeds," "))\n")
    write(outSNP, "id breed $(join(g.Breeds," "))\n")

    for tempInd in keys(g.Pop)
        yHap = zeros(Float64, g.totalNumberCommon)
        ySNP = zeros(Float64, g.numberSNPs)
        for j = 1:length(g.Segments)
            begIdx = g.SNPInfoMap[g.Segments[j].begSNP].index
            endIdx = g.SNPInfoMap[g.Segments[j].endSNP].index
            tempHap = g.Pop[tempInd].matHap[begIdx:endIdx]
            if in(tempHap, g.Segments[j].commonHap)
                rowIdx = g.Segments[j].uniqueHap[tempHap].index
                yHap[rowIdx] += 0.5
            end
            tempHap = g.Pop[tempInd].patHap[begIdx:endIdx]
            if in(tempHap, g.Segments[j].commonHap)
                rowIdx = g.Segments[j].uniqueHap[tempHap].index
                yHap[rowIdx] += 0.5
            end
            
            for k = begIdx:endIdx
                ySNP[k] = parse(Int64, g.Pop[tempInd].matHap[k]) + parse(Int64, g.Pop[tempInd].patHap[k])
            end
        end
        g.Pop[tempInd].compHap = inv(XHap'XHap) * XHap' * yHap
        write(outHap, "$tempInd $(g.Breeds[g.Pop[tempInd].breed]) $(join(g.Pop[tempInd].compHap," "))\n")

        ySNP /= 2.0
        g.Pop[tempInd].compSNP = inv(XSNP'XSNP) * XSNP' * ySNP
        write(outSNP, "$tempInd $(g.Breeds[g.Pop[tempInd].breed]) $(join(g.Pop[tempInd].compSNP," "))\n")
    end
    close(outHap)
    close(outSNP)
    println("Breed composition written to $outFile.compHap and $outFile.compSNP")
end

end


#test

#using HaplotypeModule

#cd("~/data/")

#myGenome = HaplotypeModule.makeGenome("AAN,CHA")
#HaplotypeModule.getBreedComposition!(myGenome, 500000, 0.1, "SNPInfo.AAN", "test")

