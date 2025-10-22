function output = read_json(fpath)
% This function read from a json file, parse its info and output an
% organized struct
% Example:
%     Input .json file content:
%          {
%               "usr": "admin",
%               "pwd": "12345678"
%          }
%     Output: 1x1 struct with
%          Field        Value
%          -------------------------
%          usr          'admin'
%          pwd          '12345678'
% Input:
%     fpath: .json file path, string
% Output:
%     output: struct of .json file content

fid = fopen(fpath);
raw = fread(fid,inf);
info = char(raw');
fclose(fid);
output = jsondecode(info);
