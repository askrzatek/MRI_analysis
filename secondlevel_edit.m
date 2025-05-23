%% Init
clear
clc

%% contrasts definition for each model
%% Contrast : definition case 1 & 2 for both REMINARY & PARKGAME

S1REAL_R                           = [1 0 0 0];
S1IMAGINARY_R                      = [0 1 0 0];
S2REAL_R                           = [0 0 1 0];
S2IMAGINARY_R                      = [0 0 0 1];

S1REAL_L                           = [1 0 0 0];
S1IMAGINARY_L                      = [0 1 0 0];
S2REAL_L                           = [0 0 1 0];
S2IMAGINARY_L                      = [0 0 0 1];

%% T contrast
% ANOVA2x2_SESSION - SESSION x TASK (RIGHT)

contrast_T.names = {
    
%     'SESSION1-SESSION2_R'
%     'SESSION2-SESSION1_R'
%     'IMAGINARY-REAL_R'
%     'REAL-IMAGINARY_R'
    
    'POS-SESSION-TASK-INTERACTION_R'
    'NEG-SESSION-TASK-INTERACTION_R'
    
    'S2REAL-S1REAL_R'
    'S2IMAGINARY-S1IMAGINARY_R'
    
    'S1REAL-S2REAL_R'
    'S1IMAGINARY-S2IMAGINARY_R'
    
    'S1REAL_R'
    'S1IMAGINARY_R'
    'S2REAL_R'
    'S2IMAGINARY_R'
    
    '(S1_REAL_R - S1_IMAGINARY_R) - (S2_REAL_R - S2_IMAGINARY_R)'
    '(S2_REAL_R - S2_IMAGINARY_R) - (S1_REAL_R - S1_IMAGINARY_R)'
    
    '(S1_IMAGINARY_R - S1_REAL_R) - (S2_IMAGINARY_R - S2_REAL_R)'
    '(S2_IMAGINARY_R - S2_REAL_R) - (S1_IMAGINARY_R - S1_REAL_R)'
    
    '(S1REAL_R - S1IMAGINARY_R)'
    '(S1IMAGINARY_R - S1REAL_R)'

    '(S2REAL_R - S2IMAGINARY_R)'
    '(S2IMAGINARY_R - S2REAL_R)'


}';

contrast_T.values = {
%     
% % [1 1 -1 -1]
% (S1REAL_R + S1IMAGINARY_R)       - (S2REAL_R + S2IMAGINARY_R)
% % [-1 -1 1 1]
% (S2REAL_R + S2IMAGINARY_R)       - (S1REAL_R + S1IMAGINARY_R)
% % [-1 1 -1 1]
% (S1IMAGINARY_R + S2IMAGINARY_R)  - (S1REAL_R + S2REAL_R)
% % [1 -1 1 -1]
% (S1REAL_R + S2REAL_R)            - (S1IMAGINARY_R + S2IMAGINARY_R)
%     
% [1 -1 -1 1]
(S1REAL_R + S2IMAGINARY_R)       - (S2REAL_R + S1IMAGINARY_R)
% [ -1 1 1 -1]
(S2REAL_R + S1IMAGINARY_R)       - (S1REAL_R + S2IMAGINARY_R)

% [-1 0 1 0]
S2REAL_R                         - S1REAL_R
% [ 0 -1 0 1]
S2IMAGINARY_R                    - S1IMAGINARY_R

% [1 0 -1 0]
S1REAL_R                         - S2REAL_R
% [0 1 0 -1]
S1IMAGINARY_R                    - S2IMAGINARY_R

S1REAL_R
S1IMAGINARY_R
S2REAL_R
S2IMAGINARY_R

(S1REAL_R - S1IMAGINARY_R) - (S2REAL_R - S2IMAGINARY_R)
(S2REAL_R - S2IMAGINARY_R) - (S1REAL_R - S1IMAGINARY_R)


(S1IMAGINARY_R - S1REAL_R) - (S2IMAGINARY_R - S2REAL_R)
(S2IMAGINARY_R - S2REAL_R) - (S1IMAGINARY_R - S1REAL_R)

(S1REAL_R - S1IMAGINARY_R)
(S1IMAGINARY_R - S1REAL_R)

(S2REAL_R - S2IMAGINARY_R)
(S2IMAGINARY_R - S2REAL_R)

}';


contrast_T.types = cat(1,repmat({'T'},[1 length(contrast_T.names)]));


%% F contrast

contrast_F.names = {
    
    'Avg_Condition_Effect_R'
    'Main_SESSION_Effect_R'
    'Main_TASK_Effect_R'
    'Main_SxT_INTERACTION_Effect_R'
    
}';
%%% CHECK THAT !!! EQUILIBRIUM aspect !!! should it not be 2* (S1REAL_R + S1IMAGINARY_R) ???

contrast_F.values = {
%Avg_Condition_Effect_R            = [1 1 1 1];
    S1REAL_R + S1IMAGINARY_R + S2REAL_R + S2IMAGINARY_R
%Main_SESSION_Effect_R             = [1 1 0 0];
    S1REAL_R                 + S1IMAGINARY_R
%Main_TASK_Effect_R                = [0 1 0 1];
    S1IMAGINARY_R            + S2IMAGINARY_R
%Main_SxT_INTERACTION_Effect_R     = [1 0 0 1];
    S1REAL_R                 + S2IMAGINARY_R

}';


contrast_F.types = cat(1,repmat({'F'},[1 length(contrast_F.names)]));


contrast_2x2Session.names  = [contrast_F.names  contrast_T.names];
contrast_2x2Session.values = [contrast_F.values contrast_T.values];
contrast_2x2Session.types  = [contrast_F.types  contrast_T.types];

%% ANOVA2x2_SESSION_CTL - SESSION x TASK (LEFT)

contrast_T.names = {
%     
%     'SESSION1-SESSION2_L'
%     'SESSION2-SESSION1_L'
%     'IMAGINARY-REAL_L'
%     'REAL-IMAGINARY_L'
%     
    'POS-SESSION-TASK-INTERACTION_L'
    'NEG-SESSION-TASK-INTERACTION_L'
    
    'S2REAL-S1REAL_L'
    'S2IMAGINARY-S1IMAGINARY_L'
    
    'S1REAL-S2REAL_L'
    'S1IMAGINARY-S2IMAGINARY_L'
    
    'S1REAL_L'
    'S1IMAGINARY_L'
    'S2REAL_L'
    'S2IMAGINARY_L'
    
    
}';

contrast_T.values = {
%     
%     % [1 1 -1 -1]
%     (S1REAL_L + S1IMAGINARY_L)       - (S2REAL_L + S2IMAGINARY_L)
%     % [-1 -1 1 1]
%     (S2REAL_L + S2IMAGINARY_L)       - (S1REAL_L + S1IMAGINARY_L)
%     % [-1 1 -1 1]
%     (S1IMAGINARY_L + S2IMAGINARY_L)  - (S1REAL_L + S2REAL_L)
%     % [1 -1 1 -1]
%     (S1REAL_L + S2REAL_L)            - (S1IMAGINARY_L + S2IMAGINARY_L)
%     
    % [ 1 -1 -1 1]
    (S1REAL_L + S2IMAGINARY_L)       - (S2REAL_L + S1IMAGINARY_L)
    % [ -1 1 1 -1]
    (S2REAL_L + S1IMAGINARY_L)       - (S1REAL_L + S2IMAGINARY_L)
    
    % [-1 0 1 0]
    S2REAL_L                         - S1REAL_L
    % [ 0 -1 0 1]
    S2IMAGINARY_L                    - S1IMAGINARY_L
    
    % [1 0 -1 0]
    S1REAL_L                         - S2REAL_L
    % [ 0 1 0 -1]
    S1IMAGINARY_L                    - S2IMAGINARY_L
    
    S1REAL_L
    S1IMAGINARY_L
    S2REAL_L
    S2IMAGINARY_L
    
}';


contrast_T.types = cat(1,repmat({'T'},[1 length(contrast_T.names)]));


%% F contrast

contrast_F.names = {
    
    'Avg_Condition_Effect_L'
    'Main_SESSION_Effect_L'
    'Main_TASK_Effect_L'
    'Main_SxT_INTERACTION_Effect_L'
    
}';
%%% CHECK THAT !!! EQUILIBRIUM aspect !!! should it not be 2* (S1REAL_L + S1IMAGINARY_L) ???

contrast_F.values = {
%Avg_Condition_Effect_L            = [1 1 1 1];
    S1REAL_L + S1IMAGINARY_L + S2REAL_L + S2IMAGINARY_L
%Main_SESSION_Effect_L             = [1 1 0 0];
    S1REAL_L                 + S1IMAGINARY_L
%Main_TASK_Effect_L                = [0 1 0 1];
    S1IMAGINARY_L            + S2IMAGINARY_L
%Main_SxT_INTERACTION_Effect_L     = [1 0 0 1];
    S1REAL_L                 + S2IMAGINARY_L
    
}';


contrast_F.types = cat(1,repmat({'F'},[1 length(contrast_F.names)]));


contrast_2x2Session_CTL.names  = [contrast_F.names  contrast_T.names];
contrast_2x2Session_CTL.values = [contrast_F.values contrast_T.values];
contrast_2x2Session_CTL.types  = [contrast_F.types  contrast_T.types];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Contrast : definition

LEFT_REAL_S1                            = [1 0 0 0];
LEFT_IMAGINARY_S1                       = [0 1 0 0];
RIGHT_REAL_S1                           = [0 0 1 0];
RIGHT_IMAGINARY_S1                      = [0 0 0 1];

LEFT_REAL_S2                            = [1 0 0 0];
LEFT_IMAGINARY_S2                       = [0 1 0 0];
RIGHT_REAL_S2                           = [0 0 1 0];
RIGHT_IMAGINARY_S2                      = [0 0 0 1];

%% T contrast

% ANOVA2x2_HAND - LATERALITY x TASK (S2)

contrast_T.names = {
%     
%     'RIGHT-LEFT_S2'
%     'LEFT-RIGHT_S2'
%     'IMAGINARY-REAL_S2'
%     'REAL-IMAGINARY_S2'
    
    'POS-LATERALITY-TASK-INTERACTION_S2'
    'NEG-LATERALITY-TASK-INTERACTION_S2'
    
    'RIGHT_REAL-LEFT_REAL_S2'
    'RIGHT_IMAGINARY-LEFT_IMAGINARY_S2'
    
    'LEFT_REAL-RIGHT_REAL_S2'
    'LEFT_IMAGINARY-RIGHT_IMAGINARY_S2'
    
    'LEFT_REAL_S2'
    'LEFT_IMAGINARY_S2'
    'RIGHT_REAL_S2'
    'RIGHT_IMAGINARY_S2'
    
    'RIGHT_REAL_S1 - RIGHT_IMAGINARY_S2'
    'RIGHT_IMAGINARY_S1 - RIGHT_REAL_S2'
    
}';

contrast_T.values = {
%     
%     % [1 1 -1 -1]      'LEFT-RIGHT_S2'
%     (LEFT_REAL_S2 + LEFT_IMAGINARY_S2)       - (RIGHT_REAL_S2 + RIGHT_IMAGINARY_S2)
%     % [-1 -1 1 1]      'RIGHT-LEFT_S2'
%     (RIGHT_REAL_S2 + RIGHT_IMAGINARY_S2)     - (LEFT_REAL_S2 + LEFT_IMAGINARY_S2)
%     % [-1 1 -1 1]      'IMAGINARY-REAL_S2'
%     (LEFT_IMAGINARY_S2 + RIGHT_IMAGINARY_S2)  - (LEFT_REAL_S2 + RIGHT_REAL_S2)
%     % [1 -1 1 -1]      'REAL-IMAGINARY_S2'
%     (LEFT_REAL_S2 + RIGHT_REAL_S2)            - (LEFT_IMAGINARY_S2 + RIGHT_IMAGINARY_S2)
    
    % [ 1 -1 -1 1]     'POS-LATERALITY-TASK-INTERACTION_S2'
    (LEFT_REAL_S2 + RIGHT_IMAGINARY_S2)       - (RIGHT_REAL_S2 + LEFT_IMAGINARY_S2)
    % [ -1 1 1 -1]     'NEG-LATERALITY-TASK-INTERACTION_S2'
    (RIGHT_REAL_S2 + LEFT_IMAGINARY_S2)       - (LEFT_REAL_S2 + RIGHT_IMAGINARY_S2)
    
    % [-1 0 1 0]       'RIGHT_REAL-LEFT_REAL_S2'
    RIGHT_REAL_S2                         - LEFT_REAL_S2
    % [ 0 -1 0 1]      'RIGHT_IMAGINARY-LEFT_IMAGINARY_S2'
    RIGHT_IMAGINARY_S2                    - LEFT_IMAGINARY_S2
    
    % [1 0 -1 0]       'LEFT_REAL-RIGHT_REAL_S2'
    LEFT_REAL_S2                         - RIGHT_REAL_S2
    % [ 0 1 0 -1]      'LEFT_IMAGINARY-RIGHT_IMAGINARY_S2'
    LEFT_IMAGINARY_S2                    - RIGHT_IMAGINARY_S2
    
    
    LEFT_REAL_S2
    LEFT_IMAGINARY_S2
    RIGHT_REAL_S2
    RIGHT_IMAGINARY_S2
    
     % [0 0 1 -1]
    RIGHT_REAL_S2 - RIGHT_IMAGINARY_S2
    
    % [0 0 -1 1]
    RIGHT_IMAGINARY_S2 - RIGHT_REAL_S2
    
}';


contrast_T.types = cat(1,repmat({'T'},[1 length(contrast_T.names)]));


%% F contrast

contrast_F.names = {
    
    'Avg_Condition_Effect_S2'
    'Main_LATERALITY_Effect_S2'
    'Main_TASK_Effect_S2'
    'Main_LxT_INTERACTION_Effect_S2'
    
}';
%%% CHECK THAT !!! EQUILIBRIUM aspect !!! should it not be 2* (LEFT_REAL_S2 + LEFT_IMAGINARY_S2) ???

contrast_F.values = {
    %Avg_Condition_Effect_S2            = [1 1 1 1];
    LEFT_REAL_S2 + LEFT_IMAGINARY_S2 + RIGHT_REAL_S2 + RIGHT_IMAGINARY_S2
    %Main_LATERALITY_Effect_S2             = [1 1 0 0];
    LEFT_REAL_S2                 + LEFT_IMAGINARY_S2
    %Main_TASK_Effect_S2                = [0 1 0 1];
    LEFT_IMAGINARY_S2            + RIGHT_IMAGINARY_S2
    %Main_SxT_INTERACTION_Effect_S2     = [1 0 0 1];
    LEFT_REAL_S2                 + RIGHT_IMAGINARY_S2
    
}';




contrast_F.types = cat(1,repmat({'F'},[1 length(contrast_F.names)]));


contrast_2x2Hand.names  = [contrast_F.names  contrast_T.names];
contrast_2x2Hand.values = [contrast_F.values contrast_T.values];
contrast_2x2Hand.types  = [contrast_F.types  contrast_T.types];


%% ANOVA2x2_HAND_CTL : LATERALITY x TASK (S1)

contrast_T.names = {
%     
%     'LEFT-RIGHT_S1'
%     'RIGHT-LEFT_S1'
%     'IMAGINARY-REAL_S1'
%     'REAL-IMAGINARY_S1'
    
    'POS-LATERALITY-TASK-INTERACTION_S1'
    'NEG-LATERALITY-TASK-INTERACTION_S1'
    
    'RIGHT_REAL-LEFT_REAL_S1'
    'RIGHT_IMAGINARY-LEFT_IMAGINARY_S1'
    
    'LEFT_REAL-RIGHT_REAL_S1'
    'LEFT_IMAGINARY-RIGHT_IMAGINARY_S1'
    
    
    'LEFT_REAL_S1'
    'LEFT_IMAGINARY_S1'
    'RIGHT_REAL_S1'
    'RIGHT_IMAGINARY_S1'
    
    'RIGHT_REAL_S1 - RIGHT_IMAGINARY_S1'
    'RIGHT_IMAGINARY_S1 - RIGHT_REAL_S1'
    
}';

contrast_T.values = {
%     
%     % [1 1 -1 -1]      'LEFT-RIGHT_S1'
%     (LEFT_REAL_S1 + LEFT_IMAGINARY_S1)       - (RIGHT_REAL_S1 + RIGHT_IMAGINARY_S1)
%     % [-1 -1 1 1]      'RIGHT-LEFT_S1'
%     (RIGHT_REAL_S1 + RIGHT_IMAGINARY_S1)       - (LEFT_REAL_S1 + LEFT_IMAGINARY_S1)
%     % [-1 1 -1 1]      'IMAGINARY-REAL_S1'
%     (LEFT_IMAGINARY_S1 + RIGHT_IMAGINARY_S1)  - (LEFT_REAL_S1 + RIGHT_REAL_S1)
%     % [1 -1 1 -1]      'REAL-IMAGINARY_S1'
%     (LEFT_REAL_S1 + RIGHT_REAL_S1)            - (LEFT_IMAGINARY_S1 + RIGHT_IMAGINARY_S1)
    
    % [ 1 -1 -1 1]     'POS-LATERALITY-TASK-INTERACTION_S1'
    (LEFT_REAL_S1 + RIGHT_IMAGINARY_S1)       - (RIGHT_REAL_S1 + LEFT_IMAGINARY_S1)
    % [ -1 1 1 -1]     'NEG-LATERALITY-TASK-INTERACTION_S1'
    (RIGHT_REAL_S1 + LEFT_IMAGINARY_S1)       - (LEFT_REAL_S1 + RIGHT_IMAGINARY_S1)
    
    % [-1 0 1 0]       'RIGHT_REAL-LEFT_REAL_S1'
    RIGHT_REAL_S1                         - LEFT_REAL_S1
    % [ 0 -1 0 1]      'RIGHT_IMAGINARY-LEFT_IMAGINARY_S1'
    RIGHT_IMAGINARY_S1                    - LEFT_IMAGINARY_S1
    
    % [1 0 -1 0]       'LEFT_REAL-RIGHT_REAL_S1'
    LEFT_REAL_S1                         - RIGHT_REAL_S1
    % [ 0 1 0 -1]      'LEFT_IMAGINARY-RIGHT_IMAGINARY_S1'
    LEFT_IMAGINARY_S1                    - RIGHT_IMAGINARY_S1
    
    LEFT_REAL_S1
    LEFT_IMAGINARY_S1
    RIGHT_REAL_S1
    RIGHT_IMAGINARY_S1
    
    % [0 0 1 -1]
    RIGHT_REAL_S1 - RIGHT_IMAGINARY_S1
    
    % [0 0 -1 1]
    RIGHT_IMAGINARY_S1 - RIGHT_REAL_S1
    
}';


contrast_T.types = cat(1,repmat({'T'},[1 length(contrast_T.names)]));


%% F contrast

contrast_F.names = {
    
    'Avg_Condition_Effect_S1'
    'Main_LATERALITY_Effect_S1'
    'Main_TASK_Effect_S1'
    'Main_SxT_INTERACTION_Effect_S1'
    
}';
%%% CHECK THAT !!! EQUILIBRIUM aspect !!! should it not be 2* (LEFT_REAL_S1 + LEFT_IMAGINARY_S1) ???

contrast_F.values = {
    %Avg_Condition_Effect_S1            = [1 1 1 1];
    LEFT_REAL_S1 + LEFT_IMAGINARY_S1 + RIGHT_REAL_S1 + RIGHT_IMAGINARY_S1
    %Main_LATERALITY_Effect_S1             = [1 1 0 0];
    LEFT_REAL_S1                 + LEFT_IMAGINARY_S1
    %Main_TASK_Effect_S1                = [0 1 0 1];
    LEFT_IMAGINARY_S1            + RIGHT_IMAGINARY_S1
    %Main_SxT_INTERACTION_Effect_S1     = [1 0 0 1];
    LEFT_REAL_S1                 + RIGHT_IMAGINARY_S1
    
}';


contrast_F.types = cat(1,repmat({'F'},[1 length(contrast_F.names)]));


contrast_2x2Hand_CTL.names  = [contrast_F.names  contrast_T.names];
contrast_2x2Hand_CTL.values = [contrast_F.values contrast_T.values];
contrast_2x2Hand_CTL.types  = [contrast_F.types  contrast_T.types];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% T contrast
% ANOVA2x2_TASK_IMA - LATERALITY x SESSION (IMAGINARY)
% definition contrasts

LEFT_S1_IMA                       = [1 0 0 0];
LEFT_S2_IMA                       = [0 1 0 0];
RIGHT_S1_IMA                      = [0 0 1 0];
RIGHT_S2_IMA                      = [0 0 0 1];

contrast_T.names = {
% 
%     
%     'LEFT-RIGHT_IMA'
%     'RIGHT-LEFT_IMA'
%     'SESSION1-SESSION2_IMA'
%     'SESSION2-SESSION1_IMA'

    
    'POS-HAND-SESSION-INTERACTION_IMA'
    'NEG-HAND-SESSION-INTERACTION_IMA'

    
    'RIGHT_S1-LEFT_S1_IMA'
    'RIGHT_S2-LEFT_S2_IMA'

    
    'LEFT_S1-RIGHT_S1_IMA'
    'LEFT_S2-RIGHT_S2_IMA'
    
    'LEFT_S1_IMA'
    'LEFT_S2_IMA'
    'RIGHT_S1_IMA'
    'RIGHT_S2_IMA'

    
}';

contrast_T.values = {
% 
%     
%     % [1 1 -1 -1]
%     (LEFT_S1_IMA + LEFT_S2_IMA)       - (RIGHT_S1_IMA + RIGHT_S2_IMA)
%     % [-1 -1 1 1]
%     (RIGHT_S1_IMA + RIGHT_S2_IMA)     - (LEFT_S1_IMA + LEFT_S2_IMA)
%     % [1 -1 1 -1]
%     (LEFT_S1_IMA + RIGHT_S1_IMA)      - (LEFT_S2_IMA + RIGHT_S2_IMA)
%     % [-1 1 -1 1]
%     (LEFT_S2_IMA + RIGHT_S2_IMA)      - (LEFT_S1_IMA + RIGHT_S1_IMA)

    

    
    % [1 -1 -1 1]
    (LEFT_S1_IMA + RIGHT_S2_IMA)      - (LEFT_S2_IMA + RIGHT_S1_IMA)
    % [ -1 1 1 -1]
    (LEFT_S2_IMA + RIGHT_S1_IMA)      - (LEFT_S1_IMA + RIGHT_S2_IMA)

    
    % [-1 0 1 0]
    RIGHT_S1_IMA                      - LEFT_S1_IMA
    % [ 0 -1 0 1]
    RIGHT_S2_IMA                      - LEFT_S2_IMA

    
    % [1 0 -1 0]
    LEFT_S1_IMA                       - RIGHT_S1_IMA
    % [0 1 0 -1]
    LEFT_S2_IMA                       - RIGHT_S2_IMA

    
    LEFT_S1_IMA
    LEFT_S2_IMA
    RIGHT_S1_IMA
    RIGHT_S2_IMA

    
}';


contrast_T.types = cat(1,repmat({'T'},[1 length(contrast_T.names)]));


%% F contrast

contrast_F.names = {

    
    'Avg_Condition_Effect_IMA'
    'Main_SESSION_Effect_IMA'
    'Main_HAND_Effect_IMA'
    'Main_HxS_INTERACTION_Effect_IMA'

    
}';
%%% CHECK THAT !!! EQUILIBRIUM aspect !!! should it not be 2* (S1REAL_R + S1IMAGINARY_R) ???

contrast_F.values = {
    %Avg_Condition_Effect_R                     = [1 1 1 1];
    LEFT_S1_IMA + LEFT_S2_IMA + RIGHT_S1_IMA + RIGHT_S2_IMA
    %Main_HAND_Effect_IMA                       = [1 1 0 0];
    LEFT_S1_IMA                              + LEFT_S2_IMA
    %Main_SESSION_Effect_IMA                    = [0 1 0 1];
    LEFT_S2_IMA                              + RIGHT_S2_IMA
    %Main_HxS_INTERACTION_Effect_IMA              = [1 0 0 1];
    LEFT_S1_IMA                              + RIGHT_S2_IMA

    
}';


contrast_F.types = cat(1,repmat({'F'},[1 length(contrast_F.names)]));


contrast_2x2_TASK_IMA.names  = [contrast_F.names  contrast_T.names];
contrast_2x2_TASK_IMA.values = [contrast_F.values contrast_T.values];
contrast_2x2_TASK_IMA.types  = [contrast_F.types  contrast_T.types];

%% T contrast
% ANOVA2x2_TASK_REAL - LATERALITY x SESSION (REAL)
% definition contrasts

LEFT_S1_REAL                       = [1 0 0 0];
LEFT_S2_REAL                       = [0 1 0 0];
RIGHT_S1_REAL                      = [0 0 1 0];
RIGHT_S2_REAL                      = [0 0 0 1];

contrast_T.names = {

%     
%     'LEFT-RIGHT_REAL'
%     'RIGHT-LEFT_REAL'
%     'SESSION1-SESSION2_REAL'
%     'SESSION2-SESSION1_REAL'

    
    'POS-HAND-SESSION-INTERACTION_REAL'
    'NEG-HAND-SESSION-INTERACTION_REAL'

    
    'RIGHT_S1-LEFT_S1_REAL'
    'RIGHT_S2-LEFT_S2_REAL'

    
    'LEFT_S1-RIGHT_S1_REAL'
    'LEFT_S2-RIGHT_S2_REAL'
    
    'LEFT_S1_REAL'
    'LEFT_S2_REAL'
    'RIGHT_S1_REAL'
    'RIGHT_S2_REAL'

    
}';

contrast_T.values = {

%     
%     % [1 1 -1 -1]
%     (LEFT_S1_REAL + LEFT_S2_REAL)       - (RIGHT_S1_REAL + RIGHT_S2_REAL)
%     % [-1 -1 1 1]
%     (RIGHT_S1_REAL + RIGHT_S2_REAL)     - (LEFT_S1_REAL + LEFT_S2_REAL)
%     % [1 -1 1 -1]
%     (LEFT_S1_REAL + RIGHT_S1_REAL)      - (LEFT_S2_REAL + RIGHT_S2_REAL)
%     % [-1 1 -1 1]
%     (LEFT_S2_REAL + RIGHT_S2_REAL)      - (LEFT_S1_REAL + RIGHT_S1_REAL)

    

    
    % [1 -1 -1 1]
    (LEFT_S1_REAL + RIGHT_S2_REAL)      - (LEFT_S2_REAL + RIGHT_S1_REAL)
    % [ -1 1 1 -1]
    (LEFT_S2_REAL + RIGHT_S1_REAL)      - (LEFT_S1_REAL + RIGHT_S2_REAL)

    
    % [-1 0 1 0]
    RIGHT_S1_REAL                      - LEFT_S1_REAL
    % [ 0 -1 0 1]
    RIGHT_S2_REAL                     - LEFT_S2_REAL

    
    % [1 0 -1 0]
    LEFT_S1_REAL                       - RIGHT_S1_REAL
    % [0 1 0 -1]
    LEFT_S2_REAL                       - RIGHT_S2_REAL

    
    LEFT_S1_REAL
    LEFT_S2_REAL
    RIGHT_S1_REAL
    RIGHT_S2_REAL
    
    
}';


contrast_T.types = cat(1,repmat({'T'},[1 length(contrast_T.names)]));


%% F contrast

contrast_F.names = {

    
    'Avg_Condition_Effect_REAL'
    'Main_SESSION_Effect_REAL'
    'Main_HAND_Effect_REAL'
    'Main_HxS_INTERACTION_Effect_REAL'

    
}';
%%% CHECK THAT !!! EQUILIBRIUM aspect !!! should it not be 2* (S1REAL_R + S1IMAGINARY_R) ???

contrast_F.values = {
    %Avg_Condition_Effect_R                     = [1 1 1 1];
    LEFT_S1_REAL + LEFT_S2_REAL + RIGHT_S1_REAL + RIGHT_S2_REAL
    %Main_HAND_Effect_REAL                       = [1 1 0 0];
    LEFT_S1_REAL                              + LEFT_S2_REAL
    %Main_SESSION_Effect_REAL                    = [0 1 0 1];
    LEFT_S2_REAL                              + RIGHT_S2_REAL
    %Main_HxS_INTERACTION_Effect_REAL              = [1 0 0 1];
    LEFT_S1_REAL                              + RIGHT_S2_REAL

    
}';


contrast_F.types = cat(1,repmat({'F'},[1 length(contrast_F.names)]));


contrast_2x2_TASK_REAL.names  = [contrast_F.names  contrast_T.names];
contrast_2x2_TASK_REAL.values = [contrast_F.values contrast_T.values];
contrast_2x2_TASK_REAL.types  = [contrast_F.types  contrast_T.types];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% save the model_names to create directories and their corresponding contrast structures

model_name     = {'ANOVA2x2_SxT_R',   'ANOVA2x2_SxT_L',   'ANOVA2x2_LxT_S2',  'ANOVA2x2_LxT_S1',   'ANOVA2x2_LxS_IMA', 'ANOVA2x2_LxS_REAL'};
model_contrast = { contrast_2x2Session,  contrast_2x2Session_CTL, contrast_2x2Hand,  contrast_2x2Hand_CTL,  contrast_2x2_TASK_IMA,       contrast_2x2_TASK_REAL};
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% fetch dirs for SPM
%% Load files from multiple folders
main_dir = fullfile('/network/iss/cenir/analyse/irm/users/anna.skrzatek','nifti_test');

% e_PARKGAME = exam(main_dir,'PARKGAME');
% e_REMINARY = exam(main_dir,'REMINARY_\w{2}_');
%e_REMINARYS1 = exam(main_dir,'REMINARY_\w{2}_.*1$');
%e_REMINARYS2 = exam(main_dir,'REMINARY_\w{2}_.*2$');
e_PARKGAMES2_a = exam(main_dir,'PARKGAME.*2_[ac]$'); % taking into account the expgroup "_a" for active and "_c" for control
e_PARKGAMES1_a = exam(main_dir,'PARKGAME.*1_[ac]$'); % taking into account the expgroup "_a" for active and "_c" for control
e_PARKGAMES1_a = e_PARKGAMES1_a(1:length(e_PARKGAMES2_a)); % waiting for all subjects to get a pair session

e_PARKGAMES2_c = exam(main_dir,'PARKGAME.*2_c$'); % taking into account the expgroup "_a" for active and "_c" for control
e_PARKGAMES1_c = exam(main_dir,'PARKGAME.*1_c$'); % taking into account the expgroup "_a" for active and "_c" for control


e = {e_PARKGAMES1_a, e_PARKGAMES2_a, e_PARKGAMES1_c, e_PARKGAMES2_c};
dirstat = r_mkdir(main_dir, 'secondlevel_sts_tapas_PARK');
%dirgroup = r_mkdir(char(dirstat), {'PARKGAME_a', 'PARKGAME_c'});

dirgroup = r_mkdir(char(dirstat), {'PARKGAME_all', 'PARKGAME_ac'});

done = 0;
%done = job_con_smooth('s',4); % comment this line if you don't want to smooth your contrast data
if done ==1
    for i = 1:length(e)
        %e{i}.explore
        %'REMINARY_\w{2}_.*1$'
        e{i}.addSerie('smodel_ts_tapas$','contrasts',1)

        e{i}.getSerie('contrasts').addVolume('^scon_0007','REAL_L',1)
        e{i}.getSerie('contrasts').addVolume('^scon_0008','REAL_R',1)
        e{i}.getSerie('contrasts').addVolume('^scon_0009','IMA_L',1)
        e{i}.getSerie('contrasts').addVolume('^scon_0010','IMA_R',1)

        [ec_1st, ei_1st] = e{i}.removeIncomplete;
        e{i} = ec_1st;

        e{i}.explore
        %e{2}.explore
       dirfig = 'smooth_figures_2ndlvl_ts_tapas';
    end
else
    for i = 1:length(e)
        %e{i}.explore
        %'REMINARY_\w{2}_.*1$'
        e{i}.addSerie('smodel_ts_tapas$','contrasts',1)

        e{i}.getSerie('contrasts').addVolume('^con_0008','REAL_L',1)
        e{i}.getSerie('contrasts').addVolume('^con_0009','REAL_R',1)
        e{i}.getSerie('contrasts').addVolume('^con_0010','IMA_L',1)
        e{i}.getSerie('contrasts').addVolume('^con_0011','IMA_R',1)

        [ec_1st, ei_1st] = e{i}.removeIncomplete;
        e{i} = ec_1st;

        e{i}.explore
        %e{2}.explore
       dirfig = 'figures_2ndlvl_ts_tapas';
    end
end


par.fake = 0;
par.redo = 0;
par.verbose = 2;

%% create a result directories
%dirout = r_mkdir(dirstat{1}, model_name);

%% Probably a MODEL LOOP will start here with model_name and model_contrast varying with iterations
    %% Fetch onset
    % before adding a SPM.mat to the exam we need to create one with a batch for each model with file matrices we would have fetched
    % --> modify the script job_first_level_specify(dir_func, model_dir,par) where we can fetch all needed scans according to each model, create respective folders
for group=1:2 %:4 % if only one group ready for analysis then the 4 unnecessary
    %for session= 1:2:4
    for imod=1:length(model_name)
        if group == 1 % keep if we still need 2 exams per group
            dirout = r_mkdir(dirgroup{group}, model_name);
        else
            dirout = r_mkdir(dirgroup{2}, model_name); % works only if we have no more than 2 groups ergo 2 protocols
        end
        model_dir = cellstr(dirout{imod});
        switch imod                                   
            case 1
                fact1 = 'Session';
                fact2 = 'Task';
                l11 = e{group}.getSerie('contrasts').getVolume('REAL_R').toJob;
                l12 = e{group}.getSerie('contrasts').getVolume('IMA_R').toJob;
                l21 = e{group+1}.getSerie('contrasts').getVolume('REAL_R').toJob;
                l22 = e{group+1}.getSerie('contrasts').getVolume('IMA_R').toJob;
            case 2
                fact1 = 'Session';
                fact2 = 'Task';
                l11 = e{group}.getSerie('contrasts').getVolume('REAL_L').toJob;
                l12 = e{group}.getSerie('contrasts').getVolume('IMA_L').toJob;
                l21 = e{group+1}.getSerie('contrasts').getVolume('REAL_L').toJob;
                l22 = e{group+1}.getSerie('contrasts').getVolume('IMA_L').toJob;
            case 3
                fact1 = 'Hand';
                fact2 = 'Task';
                l11 = e{group+1}.getSerie('contrasts').getVolume('REAL_L').toJob;
                l12 = e{group+1}.getSerie('contrasts').getVolume('IMA_L').toJob;
                l21 = e{group+1}.getSerie('contrasts').getVolume('REAL_R').toJob;
                l22 = e{group+1}.getSerie('contrasts').getVolume('IMA_R').toJob;
            case 4
                fact1 = 'Hand';
                fact2 = 'Task';
                l11 = e{group}.getSerie('contrasts').getVolume('REAL_L').toJob;
                l12 = e{group}.getSerie('contrasts').getVolume('IMA_L').toJob;
                l21 = e{group}.getSerie('contrasts').getVolume('REAL_R').toJob;
                l22 = e{group}.getSerie('contrasts').getVolume('IMA_R').toJob;
            case 5
                fact1 = 'Hand';
                fact2 = 'Session';
                l11 = e{group}.getSerie('contrasts').getVolume('IMA_L').toJob;
                l12 = e{group+1}.getSerie('contrasts').getVolume('IMA_L').toJob;
                l21 = e{group}.getSerie('contrasts').getVolume('IMA_R').toJob;
                l22 = e{group+1}.getSerie('contrasts').getVolume('IMA_R').toJob;
            case 6
                fact1 = 'Hand';
                fact2 = 'Session';
                l11 = e{group}.getSerie('contrasts').getVolume('REAL_L').toJob;
                l12 = e{group+1}.getSerie('contrasts').getVolume('REAL_L').toJob;
                l21 = e{group}.getSerie('contrasts').getVolume('REAL_R').toJob;
                l22 = e{group+1}.getSerie('contrasts').getVolume('REAL_R').toJob;
        end 
        addpath '/network/iss/cenir/analyse/irm/users/anna.skrzatek/'
        %% Job define model

            par.fake = 0;
            par.redo = 0;
            par.verbose = 2;
            par.run = 1;
%            par.file_reg = '^con'; % it doesn't mean anything - you can put whatever you want, it's l11, l12 etc that define the files
%            par.mask_thr = 0.07;
            par.mask_thr = 0.1;

            job_second_level_specify(fact1,fact2,l11,l12,l21,l22,model_dir,par);

    %% Estimate
            fspm = fullfile(model_dir,'SPM.mat');
            job_second_level_estimate(fspm);

    %% Contrast : write

            par.run = 1;
            par.display = 0;

            % par.sessrep = 'both';
            par.sessrep = 'none';

            par.delete_previous = 1;
            %par.report          = 1; % error on non designating spm mode
            %(?)
            par.report          = 0;
            
            job_second_level_contrast(fspm,model_contrast{imod},par);
    %% create a folder for figures before creating figures with MRIcroGL
            mkdir(dirout{imod},dirfig);
    end
end
%save e

%% Display

%!linux command for mricrogl script
%!/network/iss/cenir/software/irm/mricrogl_lx/MRIcroGL '/home/anna.skrzatek/data/nifti_test/visual_second_level_PARK.gls'