function [Header]=ReadMCSTiffInfo_2P(fname)
%extract header information from a MCS file saved as a *tiff using
%mpview tiff exporter
%PJD 5/2017
rawinfo=imfinfo(fname);
the_strings=textscan(rawinfo(1).ImageDescription,'%s','Delimiter', ':');
Header.Filename=rawinfo(1).Filename;
Header.Frame_Count=length(rawinfo);
Header.Scan_Mode=the_strings{1}{17};%
Header.X_Position = the_strings{1}{39};%
Header.Y_Position = the_strings{1}{41};%
Header.Z_Position = the_strings{1}{43};%
Header.Stack_Count = the_strings{1}{29};%
Header.Z_Interval = the_strings{1}{27};%
Header.Averaging_Count = the_strings{1}{31};%
%Header.Repeat_Count = the_strings{1}{20};
Header.Magnification = the_strings{1}{37};%
Header.Rotation = the_strings{1}{35};%
Header.X_Frame_Offset = the_strings{1}{23};%
Header.Y_Frame_Offset = the_strings{1}{25};%
Header.Frame_Duration = the_strings{1}{45};%

Header.Frame_Rate =1/str2num(Header.Frame_Duration(1:end-2));
Header.Filename = rawinfo(1).Filename;
Header.Frame_Width = rawinfo(1).Width;
Header.Frame_Height = rawinfo(1).Height;
Header.Pixel_Clock= round(1/(rawinfo(1).Width*(5/4)*rawinfo(1).Height*(Header.Frame_Rate)*(.05*1e-6))); %not used
Header.Channel_Name1 = 'null';
Header.Channel_Name2 = 'null';
Header.Channel_Name3 = 'null';
Header.Channel_Name4 = 'null';
Header.Enabled1 = 'null';
Header.Enabled2 = 'null';
Header.Enabled3 = 'null';
Header.Enabled4 = 'null';
Header.Input_Range1 = 'null';
Header.Input_Range2 = 'null';
Header.Input_Range3 = 'null';
Header.Input_Range4 = 'null';
Header.Channel_Unit3 = 'null';
Header.Channel_Unit4 = 'null';
Header.Channel_Prefix3 = 'null';
Header.Channel_Prefix4 = 'null';
Header.Conversion_Factor3 = 'null';
Header.Conversion_Factor4 = 'null';
Header.Offset3 = 'null';
Header.Offset4 = 'null';
Header.Data_Point_Per_Frame3 = 'null';
Header.Data_Point_Per_Frame4 = 'null';
end