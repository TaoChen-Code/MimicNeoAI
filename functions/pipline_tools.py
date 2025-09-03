# coding=utf-8
import logging
import subprocess
import os
from datetime import datetime
import yaml

class tools:
    def __init__(self, sys_path, log_type, log_lock):
        """
        Initialize the tools class with shared resources
        
        Args:
            sys_path (str): Base system path for file operations
            log_type (str): Type identifier for logging
            log_lock (threading.Lock): Lock for thread-safe logging
        """
        # Initialize status tracking dictionaries
        self.failed_cmds_allSamples = {}
        self.done_cmds_allSamples = {}
        self.run_cmds_allSamples = {}
        self.samples = []
        self.sys_path = sys_path + '/'
        self.log_type = log_type
        self.log_lock = log_lock
        
        # Configure logging paths and initialize logging system
        self.start_date = datetime.now().strftime("%Y-%m-%d_%H:%M:%S")
        self.start_log = self.sys_path + f'/log/{self.start_date}/{self.log_type}_start_{self.start_date}.log'
        self.cmd_log = self.sys_path + f'/log/{self.start_date}/detail/'
        
        # Create directory structure
        self.mkdir(self.sys_path+f'/log/{self.start_date}')
        
        # Configure logging handlers
        self.logger = logging.getLogger('my_logger')
        self.logger.setLevel(logging.DEBUG)
        
        # Console handler configuration
        console_handler = logging.StreamHandler()
        console_handler.setLevel(logging.DEBUG)
        
        # File handler configuration
        file_handler = logging.FileHandler(self.start_log)
        file_handler.setLevel(logging.DEBUG)
        
        # Unified formatter for both handlers
        formatter = logging.Formatter('%(asctime)s - %(levelname)s - %(message)s')
        console_handler.setFormatter(formatter)
        file_handler.setFormatter(formatter)
        
        # Add handlers to logger
        self.logger.addHandler(console_handler)
        self.logger.addHandler(file_handler)
    
    def exec_cmd(self, cmd, sample, flag = None):
        """
        Execute a system command and track its execution status
        
        Args:
            cmd (str): Command string to execute
            sample (str): Sample identifier for logging
        """
        # Parse command name for logging
        if "/" in cmd.split(" ")[0]:
            cmd_name = ''.join(cmd.split(" ")[0].split("/")[-1])
        else:
            cmd_name = ' '.join(cmd.split(" ")[0:1])

        if flag != None:
            cmd_name = cmd_name + '_' + str(flag)
            
        # Configure log file path
        current_date = datetime.now().strftime("%Y-%m-%d_%H:%M:%S")
        self.cmd_log = self.sys_path + f'/log/{self.start_date}/detail/{sample}/{current_date}_{self.log_type}_{cmd_name}.log'
        self.mkdir(self.sys_path + f'/log/{self.start_date}/detail/{sample}/')
        
        # Update running status
        with self.log_lock:
            self.status(self.run_cmds_allSamples, sample, f"{cmd_name}", 'run', '', cmd)
            
        # Execute command and track runtime
        start = datetime.now()
        with open(self.cmd_log, "a") as logfile:
            result = subprocess.call(cmd, shell=True, stdout=logfile, stderr=logfile)
        end = datetime.now()
        
        # Update status after completion
        self.run_cmds_allSamples[sample].remove(f"{cmd_name}")
        run_time = self.print_time(end-start)
        
        # Handle command result
        if result:
            self.status(self.failed_cmds_allSamples, sample, f"{cmd_name}", 'failed', run_time, cmd)
        else:
            self.status(self.done_cmds_allSamples, sample, f"{cmd_name}", 'done', run_time, cmd)

    def judge_then_exec(self, sample, cmd, file, flag = None):
        """
        Conditionally execute command based on file existence and content
        
        Args:
            sample (str): Sample identifier for logging
            cmd (str): Command to execute if conditions met
            file (str): File path to check for existence/content
        """
        if not os.path.exists(file):
            self.exec_cmd(cmd, sample, flag)
        elif os.path.isfile(file) and os.path.getsize(file) == 0:
            self.exec_cmd(cmd, sample, flag)
        else:
            self.logger.warn(file+' already exists!')

    def exec_cmd_with_time(self, cmd, sample, timeout=3600):
        """
        Execute command with timeout protection
        
        Args:
            cmd (str): Command string to execute
            sample (str): Sample identifier for logging
            timeout (int, optional): Maximum execution time in seconds
        """
        # Similar setup to exec_cmd but with timeout handling
        if "/" in cmd.split(" ")[0]:
            cmd_name = ''.join(cmd.split(" ")[0].split("/")[-1])
        else:
            cmd_name = ' '.join(cmd.split(" ")[0:1])
            
        current_date = datetime.now().strftime("%Y-%m-%d_%H:%M:%S")
        self.cmd_log = self.sys_path + f'/log/{self.start_date}/detail/{sample}/{current_date}_{self.log_type}_{cmd_name}.log'
        self.mkdir(self.sys_path + f'/log/{self.start_date}/detail/{sample}/')
        
        self.status(self.run_cmds_allSamples, sample, f"{cmd_name} timeout={timeout}", 'run', '', cmd)
        start = datetime.now()
        
        try:
            with open(self.cmd_log, "a") as logfile:
                result = subprocess.run(cmd, shell=True, stdout=logfile, stderr=logfile, timeout=timeout)
            end = datetime.now()
            run_time = self.print_time(end - start)
            
            if result.returncode != 0:
                self.status(self.failed_cmds_allSamples, sample, f"{cmd_name}", 'failed', run_time, cmd)
            else:
                self.status(self.done_cmds_allSamples, sample, f"{cmd_name}", 'done', run_time, cmd)
                
        except subprocess.TimeoutExpired:
            end = datetime.now()
            run_time = self.print_time(end - start)
            self.status(self.failed_cmds_allSamples, sample, f"{cmd_name}", 'timeout', run_time, cmd)
            self.logger.error(f"Command for sample {sample} exceeded timeout of {timeout} seconds.")
            
        self.run_cmds_allSamples[sample].remove(f"{cmd_name} timeout={timeout}")

    def judge_then_exec_with_time(self, sample, cmd, file, timeout=3600):
        """
        Conditionally execute command based on file existence and content with timeout protection
        Args:
            sample (str): Sample identifier for logging
            cmd (str): Command to execute if conditions met
            file (str): File path to check for existence/content
            timeout (int, optional): Maximum execution time in seconds
        """
        if not os.path.exists(file) or (os.path.isfile(file) and os.path.getsize(file) == 0):
            self.exec_cmd_with_time(cmd, sample, timeout)
        else:
            self.logger.warn(file + ' already exists!')
    

    def sharing_variable(self, manager, samples):
        """
        Initialize shared variables for multiprocessing
        
        Args:
            manager (multiprocessing.Manager): Manager for shared objects
            samples (list): List of sample identifiers
        """
        samples = samples + ['summary']
        # Initialize manager dictionaries
        self.failed_cmds_allSamples = manager.dict()
        self.done_cmds_allSamples = manager.dict()
        self.run_cmds_allSamples = manager.dict()
        
        # Create list containers for each sample
        for sample in samples:
            self.failed_cmds_allSamples[sample] = manager.list()
            self.done_cmds_allSamples[sample] = manager.list()
            self.run_cmds_allSamples[sample] = manager.list()

    def status(self, cmds_dict, sample, info, info_type, run_time, cmd):
        """
        Update command status and log progress
        
        Args:
            cmds_dict (dict): Dictionary to update
            sample (str): Sample identifier
            info (str): Status information
            info_type (str): Status type (run/done/failed/timeout)
            run_time (str): Formatted runtime string
            cmd (str): Original command string
        """
        # Update status dictionary
        if run_time == "0h0m0s":
            cmds_dict[sample].append(info)
        else:
            if info_type == 'run':
                cmds_dict[sample].append(info)
            else:
                cmds_dict[sample].append(info+f" {run_time}")
                
        # Generate status message
        message = ''
        for key in cmds_dict.keys():
            if len(cmds_dict[key]) != 0: 
                message += f"{key}-{cmds_dict[key][:]}\n"
        message = message[:-1]
        
        # Log based on status type
        if info_type == 'run':
            self.logger.info(f"Run {sample} {info}! Progress:{self.samples.index(sample)+1}/{len(self.samples)}. \nDetail cmd:\n{cmd}. \nRunning:\n{message}. ")
        elif info_type == 'failed':
            self.logger.error(f"{sample} {info} failed! Time: {run_time}. Please check the corresponding log file!")
        elif info_type == 'timeout':
            self.logger.error(f"{sample} {info} exceeded timeout! Time: {run_time}.")
        elif info_type == 'done':
            self.logger.info(f"{sample} {info} succeeded! Time: {run_time}. ")
        else:
            self.logger.info("status others")

    def summary(self):
        """Generate final summary of all command executions"""
        self.logger.info("end")
        # Compile completion reports
        done_message = '\n'.join([f"{key}-{self.done_cmds_allSamples[key][:]}" 
                                for key in self.done_cmds_allSamples.keys() 
                                if len(self.done_cmds_allSamples[key]) != 0])
                                
        failed_message = '\n'.join([f"{key}-{self.failed_cmds_allSamples[key][:]}" 
                                  for key in self.failed_cmds_allSamples.keys() 
                                  if len(self.failed_cmds_allSamples[key]) != 0])
        
        # Log final statuses
        if done_message:
            self.logger.info(f"all done cmds:\n{done_message}")
        if failed_message:
            self.logger.error(f"all failed cmds:\n{failed_message}")

    def write_log(self, text, type_log):
        """
        Write message to log with specified level
        
        Args:
            text (str): Message to log
            type_log (str): Log level (info/error/warn)
        """
        if type_log == "info":
            self.logger.info(text)
        elif type_log == "error":
            self.logger.error(text)
        elif type_log == 'warn':
            self.logger.warn(text)

    def mkdir(self, path):
        """
        Create directory if it doesn't exist
        
        Args:
            path (str): Directory path to create
        """
        if not os.path.exists(path):
            cmd = f"mkdir -p {path}"
            os.system(cmd)

    def print_pool_error(self, value):
        """
        Handle multiprocessing pool errors
        
        Args:
            value: Error object from pool
        """
        self.logger.error(f"error: {value}")

    def cp_configure(self, file):
        """
        Backup configuration file
        
        Args:
            file: Configuration data to backup
        """
        with open(f"{self.sys_path}/log/{self.start_date}/{self.log_type}_configure_{self.start_date}.yaml","w") as f:
            yaml.dump(file,f)

    def print_time(self, time):
        """
        Format time duration for logging (supports >24h)
        
        Args:
            time (datetime.timedelta): Time duration to format
        Returns:
            str: Formatted time string (XdXhXmXs)
        """
        total_seconds = int(time.total_seconds())  # 关键：获取总秒数（包含天数）
        days = total_seconds // 86400
        remaining_seconds = total_seconds % 86400
        hour = remaining_seconds // 3600
        minutes = (remaining_seconds % 3600) // 60
        seconds = remaining_seconds % 60
        
        # 按需拼接天数（如果 days > 0）
        if days > 0:
            return f"{days}d{hour}h{minutes}m{seconds}s"
        else:
            return f"{hour}h{minutes}m{seconds}s"

    def get_configure(self, args_configure):
        """
        Load configuration from YAML file
        
        Args:
            args_configure (str): Path to configuration file
        Returns:
            dict: Loaded configuration
        """
        if args_configure is None or len(args_configure) == 0:
            with open(f"{self.sys_path}/configures/{self.log_type}_configure.yaml","r") as file:
                configure = yaml.safe_load(file)
                self.cp_configure(configure) 
        else:
            with open(f"{args_configure}","r") as file:
                configure = yaml.safe_load(file)
                self.cp_configure(configure)
        return configure

    def get_pathes(self, args_pathes):
        """
        Load pathes from YAML file
        
        Args:
            args_pathes (str): Path to pathes file
        Returns:
            dict: Loaded configuration
        """
        if args_pathes is None or len(args_pathes) == 0:
            with open(f"{self.sys_path}/configures/{self.log_type}_pathes.yaml","r") as file:
                pathes = yaml.safe_load(file)
        else:
            with open(f"{args_pathes}","r") as file:
                pathes = yaml.safe_load(file)
        return pathes

    def write2shell(self, cmd, fileName):
        """
        Write command to shell script for later execution
        
        Args:
            cmd (str): Command to save
            fileName (str): Name for the shell script
        Returns:
            str: Path to created shell script
        """
        file = f"{self.sys_path}/log/{self.start_date}/{fileName}.sh"
        with open(file,"a") as f:
            f.write(f"{cmd}\n")
        return file