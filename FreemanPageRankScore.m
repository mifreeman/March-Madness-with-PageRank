%%Author: Dr. Tim Chartier and Michael Freeman%%

function scorePageRankESPN(gameFilename, teamFilename)

%% Input parameters

%gameFilename = 'data/2015games.txt';
%teamFilename = 'data/2015teams.txt';

weighting = 2;

weightHomeWin = .7;
weightAwayWin = 1.3;
weightNeutralWin = 1;
segmentWeighting = [.5 1 1.3 1.4];

%%
selectionSundayList = {'03/10/2002','03/16/2003','03/14/2004','03/13/2005','03/12/2006',...
    '03/11/2007','03/16/2008','03/15/2009','03/14/2010','03/13/2011',...
    '03/11/2012','03/17/2013','03/16/2014','03/15/2015','03/13/2016','03/12/2017','03/11/2018'};

% Load the team names into an array
fid = fopen(teamFilename);
counter = 0;
teamname = fgetl(fid);
while (ischar(teamname))
    [token, remain] = strtok(teamname); teamname = strtok(remain);
    teamname=cellstr(teamname);
    counter = counter + 1;
    teamNames(counter) = teamname;
    teamname = fgetl(fid);
end
fclose(fid);
numTeams = counter;

%% Load the games
games=load(gameFilename);
% columns of games are:
%	column 1 = days since 1/1/0000
%	column 2 = date in YYYYMMDD format
%	column 3 = team1 index
%	column 4 = team1 homefield (1 = home, -1 = away, 0 = neutral)
%	column 5 = team1 score
%	column 6 = team2 index
%	column 7 = team2 homefield (1 = home, -1 = away, 0 = neutral)
%	column 8 = team2 score

% Pull off just the days
% days=games(:,1); days may not match
dayColumn = mod(games(:,2),100);
monthColumn = mod(floor(games(:,2)/100),100);
yearColumn = floor(games(:,2)/(100*100));
days = datenum(yearColumn,monthColumn,dayColumn);

%% Pull off year
lastDateYear = games(end,2);
yearOfData = floor(lastDateYear/10000);
selectionSunday = datenum(selectionSundayList{yearOfData - 2001});

numGames =  length(find(days <= selectionSunday));
% This is the number of games played prior to the tournament.

%% Find March Madness teams
marchMadnessTeams = [games(end,3), games(end,6)];
numTotalGames = size(games,1);
for i=numTotalGames:-1:1
    currentDay = mod(games(i,2),100);
    currentMonth = mod(floor(games(i,2)/100),100);
    currentYear = floor(games(i,2)/(100*100));
    daysAfter2000 = datenum(currentYear,currentMonth,currentDay);
    
    if (daysAfter2000 <= selectionSunday)
        break;
    end
    teamsInGame = [games(i,3),games(i,6)];
    if ~isempty(intersect(teamsInGame,marchMadnessTeams))
        marchMadnessTeams = union(teamsInGame,marchMadnessTeams);
    end
end

%% Calculate the weights

%create weight vector w 
w=zeros(numGames, 1); 

if weighting==0
    % No Weighting
    w=ones(numGames, 1); 
elseif weighting==1
    % for linear weighted time
    for i=1:numGames
        w(i)=(days(i)-days(1))/(days(numGames)-days(1)); 
    end
elseif weighting==2 % interval weighting 
    w=ones(numGames, 1); 
    for i=1:numGames
        weightIndex = ceil((days(i)-days(1)+1)/(days(numGames)-days(1)+1)*length(segmentWeighting)); 
        w(i) = segmentWeighting(weightIndex);
    end
else
   fprintf('No weighting indicated so default of uniform weighting will be assumed.\n'); 
   % No Weighting
    w=ones(numGames, 1); 
end


%% Create game network matrix

gameMatrix = zeros(numTeams);
for i=1:numGames
    team1ID = games(i, 3);
    team1Score = games(i, 5);
    team1Loc = games(i, 4);
    team2ID = games(i, 6);
    team2Score = games(i, 8);
    team2Loc = games(i, 7);
    
    if team1Score > team2Score
        % Team 1 won
        if (team1Loc == 1)       % Home win
            gameMatrix(team2ID, team1ID) = gameMatrix(team2ID, team1ID) + weightHomeWin*w(i);
        elseif (team1Loc == -1)  % Away win
            gameMatrix(team2ID, team1ID) = gameMatrix(team2ID, team1ID) + weightAwayWin*w(i);
        else                       % Neutral court win
            gameMatrix(team2ID, team1ID) = gameMatrix(team2ID, team1ID) + weightNeutralWin*w(i);
        end
    elseif team1Score < team2Score
        % Team 2 Won
        if (team2Loc == 1)       % Home win
            gameMatrix(team1ID, team2ID) = gameMatrix(team1ID, team2ID) + weightHomeWin*w(i);
        elseif (team2Loc == -1)  % Away win
            gameMatrix(team1ID, team2ID) = gameMatrix(team1ID, team2ID) + weightAwayWin*w(i);
        else                       % Neutral court win
            gameMatrix(team1ID, team2ID) = gameMatrix(team1ID, team2ID) + weightNeutralWin*w(i);
        end
    else  % it is a tie and make 1/2 a win and 1/2 a loss for both teams
        fprintf('There is a tie in the data.  Not processed.\n')
    end
end

%% Create Google matrix
alpha = 0.55; % teleportation parameter

googleMatrix = (1-alpha)/numTeams*ones(numTeams);
for i=1:numTeams
    rowSum = sum(gameMatrix(i,:));
    if (rowSum ~= 0)
        indexOfLink = find(gameMatrix(i,:) > 0);
        for j=1:length(indexOfLink)
            googleMatrix(i,indexOfLink(j)) = googleMatrix(i,indexOfLink(j)) + gameMatrix(i,indexOfLink(j))*alpha/rowSum;
        end
    else % dangling node
        %       fprintf('Dangling node = %d\n',i)
        googleMatrix(i,:) = 1/numTeams*ones(1,numTeams);
    end
end

%% Calculate ratings

G = googleMatrix^10000;
r = G(1,:);
[sortedr,index]=sort(r,'descend');

%% Find the number of correct predictions
numGamesBeforeRound1 = length(find(days <= selectionSunday+3));

predictedMadnessWins = zeros(numTeams,1);
missedMadnessWins = zeros(numTeams,1);
predictionMissed = 0;
predictionMade = 0;
for i=numGamesBeforeRound1+1:numTotalGames
    team1ID = games(i, 3);
    team1Score = games(i, 5);
    team2ID = games(i, 6);
    team2Score = games(i, 8);
    
    if (length(intersect(marchMadnessTeams,[team1ID,team2ID])) == 2)
        if (team1Score > team2Score) && (r(team1ID) > r(team2ID)) && ...
                (missedMadnessWins(team1ID) == 0)
            predictedMadnessWins(team1ID) = predictedMadnessWins(team1ID) + 1;
        elseif (team1Score > team2Score) && (r(team1ID) < r(team2ID))
            missedMadnessWins(team1ID) = missedMadnessWins(team1ID) + 1;
        elseif (team2Score > team1Score) && (r(team2ID) > r(team1ID)) && ...
                (missedMadnessWins(team2ID) == 0)
            predictedMadnessWins(team2ID) = predictedMadnessWins(team2ID) + 1;
        elseif (team2Score > team1Score) && (r(team2ID) < r(team1ID))
            missedMadnessWins(team2ID) = missedMadnessWins(team2ID) + 1;
        end
    end
end


%% Find the ESPN score
espnScore = 0;
scores = [10, 30, 70, 150, 310, 630];
teamsWithCorrectPredictions = find(predictedMadnessWins);
team_roundscore = zeros(length(teamsWithCorrectPredictions),1);
for indexOfTeam=1:length(teamsWithCorrectPredictions)
    increaseAmt = scores(predictedMadnessWins(teamsWithCorrectPredictions(indexOfTeam)));
    team_roundscore(indexOfTeam) = increaseAmt;
    espnScore = espnScore + scores(predictedMadnessWins(teamsWithCorrectPredictions(indexOfTeam)));
    %fprintf('PageRank ESPN score for %4d = %d, increased by %d\n',yearOfData,espnScore, increaseAmt);
end
%round by round scoring 
round_64 = length(team_roundscore); 
round_32 = 0; 
round_16 = 0;
round_8 = 0;
round_4 = 0; 
round_2 = 0; 

for i = 1:length(team_roundscore)
    if team_roundscore(i) >= 30 
        round_32 = round_32 + 1;
    if team_roundscore(i) >= 70 
        round_16 = round_16 + 1;
    if team_roundscore(i) >= 150 
        round_8 = round_8 + 1;
    if team_roundscore(i) >= 310 
        round_4 = round_4 + 1;
    if team_roundscore(i) >= 630
        round_2 = round_2 + 1;
    end
    end 
    end
    end
    end
end

percent_64 = round_64 / 32.0;
percent_32 = round_32 / 16.0;
percent_16 = round_16 / 8.0;
percent_8 = round_8 / 4.0;
percent_4 = round_4 / 2.0;
percent_2 = round_2 / 1.0;
    
export_m = [yearOfData; percent_64; percent_32; percent_16; percent_8; percent_4; percent_2];

dlmwrite('rounddata.csv', export_m, '-append')%writes round by round data to a csv

%% Print March Madness teams

rMarch = r(marchMadnessTeams);
teamsMarch = teamNames(marchMadnessTeams);

[sortedr,index]=sort(rMarch,'descend');

%% Statistics to Output to Screen
% fprintf('\n\n******* PageRank March Madness Rating Method *******\n\n');
% fprintf('===========================\n');
% fprintf('Results for Top %d Teams\n',length(teamsMarch));
% fprintf('===========================\n');
% fprintf('Rank   Rating    Team   \n');
% fprintf('===========================\n');
% for i=1:length(teamsMarch)
%     fprintf('%4d  %8.5f  %s\n',i,sortedr(i),cell2mat(teamsMarch(index(i))));
% %    fprintf('%4d %s\n',i,cell2mat(teamsMarch(index(i))));
% end


fprintf('\n');
%fprintf('Percent Correct in Round of 64 = %f\n', percent_64 * 100);
%fprintf('Percent Correct in Round of 32 = %f\n', percent_32 * 100);
%fprintf('Percent Correct in Round of 16 = %f\n', percent_16 * 100);
%fprintf('Percent Correct in Round of 8 = %f\n', percent_8 * 100);
%fprintf('Percent Correct in Round of 4 = %f\n', percent_4 * 100);
%fprintf('Percent Correct in Round of 2 = %f\n', percent_2 * 100);
fprintf('PageRank ESPN score for %4d = %d\n',yearOfData,espnScore);
fprintf('Percentage correct = %f\n',sum(predictedMadnessWins)/63*100)
