function [ badChannels ] = RPS_channelCheckbox()
% CARE_CHANNELCHECKBOX is a function, which displays a small GUI for the 
% selection of bad channels. It returns a vector including the numbers
% of the bad channels
%
% Use as
%   [ badChannel ] = JAI_channelCheckbox()
%
% SEE also UIFIGURE, UICHECKBOX, UIBUTTON, UIRESUME, UIWAIT

% Copyright (C) 2018, Daniel Matthes, MPI CBS

% -------------------------------------------------------------------------
% Create GUI
% -------------------------------------------------------------------------
SelectBadChannels = uifigure;
SelectBadChannels.Position = [150 400 375 215];
SelectBadChannels.Name = 'Select bad channels';

% Create Ch01CheckBox
Elec.Ch01 = uicheckbox(SelectBadChannels);
Elec.Ch01.Text = 'Ch01';
Elec.Ch01.Position = [45 150 80 15];
% Create Ch02CheckBox
Elec.Ch02 = uicheckbox(SelectBadChannels);
Elec.Ch02.Text = 'Ch02';
Elec.Ch02.Position = [125 150 80 15];
% Create Ch03CheckBox
Elec.Ch03 = uicheckbox(SelectBadChannels);
Elec.Ch03.Text = 'Ch03';
Elec.Ch03.Position = [205 150 80 15];
% Create Ch04CheckBox
Elec.Ch04 = uicheckbox(SelectBadChannels);
Elec.Ch04.Text = 'Ch04';
Elec.Ch04.Position = [285 150 80 15];

% Create Ch05CheckBox
Elec.Ch05 = uicheckbox(SelectBadChannels);
Elec.Ch05.Text = 'Ch05';
Elec.Ch05.Position = [45 125 80 15];
% Create Ch06CheckBox
Elec.Ch06 = uicheckbox(SelectBadChannels);
Elec.Ch06.Text = 'Ch06';
Elec.Ch06.Position = [125 125 80 15];
% Create Ch07CheckBox
Elec.Ch07 = uicheckbox(SelectBadChannels);
Elec.Ch07.Text = 'Ch07';
Elec.Ch07.Position = [205 125 80 15];
% Create Ch08CheckBox
Elec.Ch08 = uicheckbox(SelectBadChannels);
Elec.Ch08.Text = 'Ch08';
Elec.Ch08.Position = [285 125 80 15];

% Create Ch09CheckBox
Elec.Ch09 = uicheckbox(SelectBadChannels);
Elec.Ch09.Text = 'Ch09';
Elec.Ch09.Position = [45 100 80 15];
% Create Ch10CheckBox
Elec.Ch10 = uicheckbox(SelectBadChannels);
Elec.Ch10.Text = 'Ch10';
Elec.Ch10.Position = [125 100 80 15];
% Create Ch11CheckBox
Elec.Ch11 = uicheckbox(SelectBadChannels);
Elec.Ch11.Text = 'Ch11';
Elec.Ch11.Position = [205 100 80 15];
% Create Ch12CheckBox
Elec.Ch12 = uicheckbox(SelectBadChannels);
Elec.Ch12.Text = 'Ch12';
Elec.Ch12.Position = [285 100 80 15];

% Create Ch13CheckBox
Elec.Ch13 = uicheckbox(SelectBadChannels);
Elec.Ch13.Text = 'Ch13';
Elec.Ch13.Position = [45 75 80 15];
% Create Ch14CheckBox
Elec.Ch14 = uicheckbox(SelectBadChannels);
Elec.Ch14.Text = 'Ch14';
Elec.Ch14.Position = [125 75 80 15];
% Create Ch15CheckBox
Elec.Ch15 = uicheckbox(SelectBadChannels);
Elec.Ch15.Text = 'Ch15';
Elec.Ch15.Position = [205 75 80 15];
% Create Ch16CheckBox
Elec.Ch16 = uicheckbox(SelectBadChannels);
Elec.Ch16.Text = 'Ch16';
Elec.Ch16.Position = [285 75 80 15];

% Create SaveButton
btn = uibutton(SelectBadChannels, 'push');
btn.ButtonPushedFcn = @(btn, evt)SaveButtonPushed(SelectBadChannels);
btn.Position = [137 27 101 21];
btn.Text = 'Save';

% -------------------------------------------------------------------------
% Wait for user input and return selection after btn 'save' was pressed
% -------------------------------------------------------------------------
% Wait until btn is pushed
uiwait(SelectBadChannels);

if ishandle(SelectBadChannels)                                              % if gui still exists
  badChannels  = [ Elec.Ch01.Value; Elec.Ch02.Value; Elec.Ch03.Value; ...    % return existing selection
                  Elec.Ch04.Value; Elec.Ch05.Value; Elec.Ch06.Value; ...
                  Elec.Ch07.Value; Elec.Ch08.Value; Elec.Ch09.Value; ...
                  Elec.Ch10.Value; Elec.Ch11.Value; Elec.Ch12.Value; ...
                  Elec.Ch13.Value; Elec.Ch14.Value; Elec.Ch15.Value; ...
                  Elec.Ch16.Value ];
  number      = 1:16;
  badChannels  = number(badChannels);
  if isempty(badChannels)
    badChannels = [];
  end
  delete(SelectBadChannels);                                                % close gui
else                                                                        % if gui was already closed (i.e. by using the close symbol)
  badChannels = [];                                                          % return empty selection
end

end

% -------------------------------------------------------------------------
% Event Functions
% -------------------------------------------------------------------------
% Button pushed function: btn
function  SaveButtonPushed(SelectBadChannels)
  uiresume(SelectBadChannels);                                              % resume from wait status                                                                             
end
