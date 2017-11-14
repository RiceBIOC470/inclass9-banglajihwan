% Inclass assignment 9
%Gb comments
1) 100
2) 100
3) 100 
4) 100
Overall: 100


% The accession number for human NOTCH1 mRNA is AF308602
% 1. Read the information from this entry into matlab
genebank_dat = getgenbank ('AF308602')  
% 2. Write code that runs a blast query on the first 500 base pairs of this
% gene against the refseq_rna database
seq_begin = genebank_dat.Sequence(1:500); 
[requestID, requestTime] = blastncbi (seq_begin, 'blastn', 'Database', 'refseq_rna' ); 
blast_data = getblast(requestID, 'WaitTime', requestTime); 
% 3. Find the three highest scoring hits from other species and identify
% the length of the alignment and fraction of matches/mismatches. 

% -> I made a loop to display all the struts under HITS 
for i = 1:50
    disp(i)
    blast_data.Hits(i) 
end ; 
% -> I can see that the second HITS is already a different species.
% Therefore  I will store these struts: 
hit1 = blast_data.Hits(2) % Pan troglodytes
hit2 = blast_data.Hits(6) % Rhinopithecus biet
hit3 = blast_data.Hits(7) % Cercocebus atys

%% Pan troglodytes (PT) 
%length
length1 = hit1.HSPs.QueryIndices;
lengthPT = length1(2) - length1(1)+1; 

%ratio 
id1 = hit1.HSPs.Identities;
match1= id1.Match
unmatch1 = (id1.Possible - id1.Match) 
ratio1 = match1/unmatch1 

%% Rhinopithecus biet (Rb) 
%length
length2 = hit2.HSPs.QueryIndices;
lengthRb = length2(2) - length2(1)+1; 
%ratio
id2 = hit2.HSPs.Identities;
match2= id2.Match
unmatch2 = (id2.Possible - id2.Match) 
ratio2 = match2/unmatch2

%% Cercocebus atys (Ca) 
%length
length3 = hit3.HSPs.QueryIndices;
legnthCa = length3(2) - length3(1) +1; 
%ratio
id3 = hit3.HSPs.Identities;
match3= id3.Match
unmatch3 = (id3.Possible - id3.Match) 
ratio3 = match3/unmatch3

% 4. Run the same query against the database est_human. Comment on the
% sequences that you find. 
[requestID1, requestTime1] = blastncbi (seq_begin, 'blastn', 'Database', 'est_human' ); 
blast_data1 = getblast(requestID1, 'WaitTime', requestTime1); 
% all the sequence I get is from the human genome. 
