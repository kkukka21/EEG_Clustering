function vers = eegplugin_eeg_clustering(fig, try_strings, catch_strings)
    vers = '1.0'; % plugin version
    
    % Create menu under "Tools"
    menu = findobj(fig, 'tag', 'tools'); 
    uimenu(menu, 'label', 'EEG Clustering', ...
        'callback', 'EEG_Clustering'); % Calls your app
end
