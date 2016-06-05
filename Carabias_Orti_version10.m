[d1,sr] = audioread('k310_1_synth.wav');
[d2,sr] = audioread('k310_1_gould.wav');

d1 = d1(:,1);
d2 = d2(:,1);

nfft = 512;
window = 2048;
noverlap = window - nfft;

[s, f, t] = spectrogram (d1, window, noverlap);
[sRef, fRef, tRef] = spectrogram (d2, window, noverlap);


Y = abs(s);
Y = Y +0.000001;
X = abs(sRef);
 

%imagesc (t, f, Y);



%%
MIDIFilename = 'k310_1.mid';


% Rewrite MIDI with fixed times
nmat = readmidi_java(MIDIFilename,true);
nmat(:,7) = nmat(:,6) + nmat(:,7);
noteArray = unique(nmat(:,4));


sheetMatrix = zeros(max(noteArray) - min(noteArray), floor(nmat(length(nmat), 7) * 44100/nfft));


for i = 1 : length(nmat)
    notePitch = nmat(i,4);
    onset = floor( nmat(i,6) * 44100 / nfft) + 1;
    offset = floor ( nmat(i,7) * 44100 / nfft) +1;
    
    sheetMatrix (notePitch, onset:offset) = 1;
end


imagesc(sheetMatrix)
axis('xy')


if size(Y,2)> size(sheetMatrix,2) 
    Y(:, size(sheetMatrix,2) + 1 : size(Y,2) ) = [];
end

if size(Y,2) < size(sheetMatrix,2) 
    Y = horzcat(Y, zeros(size(Y,1), size(sheetMatrix,2) - size(Y,2))); 
end


noteCombination = zeros(10, length(sheetMatrix));

for j = 1 : length(sheetMatrix)
    noteCombinationTemp = [];
    for i = 1: size(sheetMatrix,1)
       if sheetMatrix(i,j) == 1
           noteCombinationTemp(end+1) = i;
       end   
    end
    noteCombination(1:length(noteCombinationTemp),j) =  noteCombinationTemp;
end

[C, idx] = unique(noteCombination', 'rows');

Ctime = zeros(size(noteCombination'));

for i = 1: length(C)
    Ctime(idx(i), :) = C(i,:);
end

deleteIndex = [];
for i = 2: length(Ctime)
    if Ctime(i,:) == zeros(1,11)
        deleteIndex(length(deleteIndex)+1) = i;
    end
end
%%
deleteIndex = sort(deleteIndex, 'descend');

for i = 1 : length(deleteIndex)
   Ctime(deleteIndex(i),:) = []; 
end



R = zeros(length(Ctime), length(sheetMatrix));
for i = 1 : length(sheetMatrix)
    R(find(ismember(Ctime, noteCombination(:,i)', 'rows')), i) = 1;
end

%imagesc(R)

Q = zeros(max(noteArray), length(C));

for i = 1: length(C)
    for j = 1 : size(C,2);
        if Ctime(i, j) ~= 0
            Q (Ctime(i,j), i ) = 1;
        end
    end
end

%imagesc(Q)

%% MIDI States-MIDI Notes Matrix S(n,m)

midiState = zeros(10, length(sheetMatrix));

noteCombinationOld = [0];
index = 1;

for j = 1 : length(sheetMatrix)
    noteCombinationTemp = [];
    for i = 1: size(sheetMatrix,1)
       if sheetMatrix(i,j) == 1
           noteCombinationTemp(end+1) = i;
       end   
    end
    if length(noteCombinationTemp) ~= length(noteCombinationOld)
        midiState(1:length(noteCombinationTemp),index) =  noteCombinationTemp;
        index = index + 1;
    else if noteCombinationTemp ~= noteCombinationOld 
        midiState(1:length(noteCombinationTemp),index) =  noteCombinationTemp;   
        index = index + 1;
        end
    end
    noteCombinationOld = noteCombinationTemp;
end

midiState(:, index+1:size(midiState,2)) = [];

S = zeros(max(noteArray), length(midiState));

for j = 1 : length(S)
    for i = 1:size(midiState,1)
        if midiState(i,j) ~= 0
            S(midiState(i,j), j) = 1;
        end
    end
end
midiState = midiState';
imagesc(S)    

%% States-to-combination matrix H(m,k)
H = zeros(size(midiState,1), size(Ctime,1));

[Lia,Locb] = ismember(midiState, Ctime, 'rows');


for i = 1: size(H,1)
    H(i, Locb(i)) = 1;
end        

imagesc(H)


%%
G = R;
B = rand(size(Y,1), size(G,1));
beta = 2;
Yhat = B * G + 0.000001;


%
B = B .* ((Y .* (Yhat .^(beta-2) ) * G') ./ ((Yhat .^ (beta-1)) * G'));
G = G .* ( B' * (Y .* (Yhat .^(beta-2) )) ./ (B' * (Yhat .^ (beta-1))));

B = betaNormC(B,beta);
B(find(isnan(B)))=0;
G(find(isnan(G)))=0;
%

Yhat = B * G + 0.000001;
betaDivergence = betaDivergenceMatrix(Y, Yhat, beta)

for i = 1:3
    B = B .* ((Y .* (Yhat .^(beta-2) ) * G') ./ ((Yhat .^ (beta-1)) * G'));
    B = betaNormC(B,beta);
    B(find(isnan(B)))=0;


    
    Yhat = B * G + 0.000001;

    betaDivergence = betaDivergenceMatrix(Y, Yhat, beta)
end


% Divergence Matrix Phi


PhiTest = zeros(size(B,2), size(X,2));

gainKT = betaNormC(B, beta)' .^ (beta-1) * X;

%
if beta == 1
    xTerm = repmat(sum(X),size(B,2) ,1);
    yTerm = repmat(sum(B)', 1, size(X,2)) .* gainKT;
    logTerm = (1 ./ B)' * X ./ gainKT;
    Phi = xTerm .* (log(logTerm) ./ log(10)) - xTerm + yTerm;

end



if and(beta ~= 0, beta ~= 1)
    xBetaTerm = repmat(sum(X .^ beta),size(B,2) ,1);
    yBetaTerm = repmat(sum(B .^ beta)', 1, size(X,2)) .* (gainKT .^ beta) * (beta-1);

    Phi = (xBetaTerm +  yBetaTerm - B' .^ (beta-1) * (X .* beta) .* gainKT) / beta / (beta -1);
end

%%
for t = 1: size(X,2)
    K = zeros(size(B,2),1);
    for k = 1 : 10
        g = sum (X(:,t) .* B(:,k) .^ (beta-1)) / sum(  B(:,k) .^ beta);
        K(k) = sum( ( X(:,t) .^ beta + (beta-1) * ((g * B(:,k)) .^ beta)  - beta * X(:,t) .* (g * B(:,k)) .^ (beta-1)) / ( beta * (beta -1)) );

    end
    k = find(K == min(K));

    PhiTest(:,t) = K;
end


%%

D = R' * Phi;


[p, q, DD, sc] = dpfast(D);


imagesc(D)
colormap(1-gray)
hold on; plot(q,p,'r'); hold off

%%
nmat = readmidi_java(MIDIFilename,true);
nmat(:,7) = nmat(:,6) + nmat(:,7);


% map the times
for i = 1:length(nmat)
    index = min(find(p == ceil(nmat(i,6) * (44100 / nfft))+1));
    nmat(i,6) = q(index) / 44100 * nfft;
    index2 = max(find(p == floor(nmat(i,7) * (44100 / nfft))));
    nmat(i,7) = q(index2) / 44100 * nfft;
    
end



%
% end times back to durations
nmat(:,7) = nmat(:,7) - nmat(:,6);

nmat(find(nmat(:,7) <= 0),7) = 0.08;



MFout = 'testCarabiasGould988.mid';
writemidi_seconds(nmat,MFout);

