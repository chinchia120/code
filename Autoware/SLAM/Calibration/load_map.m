function data = load_map(filename)
%     data = pcread(filename);
    data = load(filename);
    data = data.HDStripMap;
    data = data.ptCloud;
end