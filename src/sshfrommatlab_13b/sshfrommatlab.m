function channel  =  sshfrommatlab(userName,hostName,password)
%SSHFROMMATLAB connects Matlab to a remote computer via a secure shell
%
% CONN  =  SSHFROMMATLAB(USERNAME,HOSTNAME,PASSWORD)
%
% Inputs:
%   USERNAME is the user name required for the remote machine
%   HOSTNAME is the name of the remote machine
%   PASSWORD is the password for the account USERNAME@HOSTNAME
%
% Outputs:
%   CONN is a Java ch.ethz.ssh2.Connection object
%
% See also SSHFROMMATLABCLOSE, SSHFROMMATLABINSTALL, SSHFROMMATLABISSUE
%
% (c) 2008 British Oceanographic Data Centre
%    Adam Leadbetter (alead@bodc.ac.uk)
%     2010 Boston University - ECE
%    David Scott Freedman (dfreedma@bu.edu)
%    Version 1.3
%

%
%  Invocation checks
%
usesshlib = 2;
if(nargin  ~=  3)
    error('Error: SSHFROMMATLAB requires 3 input arguments...');
end
if(~ischar(userName)  || ~ischar(hostName)  ||  ~ischar(password))
    error...
        (['Error: SSHFROMMATLAB requires all input ',...
        'arguments to be strings...']);
end
%
%  Build the connection using the JSch package
%
switch usesshlib
    case 1
        try
            import ch.ethz.ssh2.*;
            try
                channel  =  Connection(hostName);
                channel.connect();
            catch
                error(['Error: SSHFROMMATLAB could not connect to the'...
                    ' remote machine %s ...'],...
                    hostName);
            end
        catch
            error('Error: SSHFROMMATLAB could not find the SSH2 java package');
        end
        %
        %  Check the authentication for login...
        %
        isAuthenticated = channel.authenticateWithPassword(userName,password);
        if(~isAuthenticated)
            error...
                (['Error: SSHFROMMATLAB could not authenticate the',...
                ' SSH connection...']);
        end
    case 2
        try
            import net.schmizz.sshj.*;
            import net.schmizz.keepalive.*;%KeepAliveProvider;
            import java.io.File;
            
            %             import net.i2p.*;
            %             import net.schmizz.sshj.DefaultConfig;
            %             import net.schmizz.sshj.SSHClient;
            import net.schmizz.sshj.connection.channel.direct.Session;
            %             import net.schmizz.sshj.transport.verification.PromiscuousVerifier;
            %             import net.schmizz.sshj.SSHClient;
            %             import net.schmizz.sshj.common.*;
            %import net.schmizz.sshj.connection.*;
            try
                channel  =  SSHClient();
                %                 khFile = File(transport.verification.OpenSSHKnownHosts.detectSSHDir(), 'known_hosts');
                %                 channel.addHostKeyVerifier(transport.verification.ConsoleKnownHostsVerifier(khFile, java.lang.System.console()));
                %                 channel.loadKnownHosts
                channel.addHostKeyVerifier(transport.verification.PromiscuousVerifier);
                
                %                 channel.authPublickey('mbeining')
                %                 channel.loadKnownHosts;
                channel.connect(hostName)
            catch
                error(['Error: SSHFROMMATLAB could not connect to the'...
                    ' remote machine %s ...'],...
                    hostName);
            end
        catch
            error('Error: SSHFROMMATLAB could not find the SSH2 java package');
        end
        %
        %  Check the authentication for login...
        
        channel.authPassword(userName,password);
        if ~(channel.isAuthenticated)
            error...
                (['Error: SSHFROMMATLAB could not authenticate the',...
                ' SSH connection...']);
        end
        
end