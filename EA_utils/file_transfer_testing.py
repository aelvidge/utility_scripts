#usage example:
#####################################################################################################
top_source_directory_path = r'C:\Users\AE78333\OneDrive - SSE PLC\Documents\Projects'
top_destination_directory_path = r'W:\01 Projects\03 Offshore'

#source_and_destination_subdirectories must be directories, not files
source_and_destination_subdirectories = [['Celtic_Sea', 'Celtic Sea']]

#objects_not_to_copy: No file or directory naned the same as any item on this list will be copied. Each item must include any file extension
objects_not_to_copy = ['Celtic_Sea_project_notes', 'Celtic Sea - Shortcut to Sharepoint project space.lnk', 'Python.lnk', 'all_ALP_output_files']
#paths_not_to_copy: List of paths (with any file extensions) to not copy across
paths_not_to_copy = []

replace_or_update_or_keep_or_forbid_existing_duplicate_files = 'u' #'r' = replace, 'u' = update (copy across if source file is newer than destination file, otherwise keep), 'k' = keep (don't copy across and continue to next file), 'f' = forbid (raise error)
#####################################################################################################
import shutil
import win32com.client # to install, run pip install pywin32
import os
from pathlib import Path
import time

st = time.time()

def recursive_selective_file_transfer(source_object_path, destination_object_path, objects_not_to_copy, paths_not_to_copy):
    source_object_contents = os.listdir(source_object_path)
    failed_actions = []
    for source_object_content in source_object_contents:
        source_object_content_path = fr'{source_object_path}\{source_object_content}'
        if source_object_content not in objects_not_to_copy and source_object_content_path not in paths_not_to_copy:
            if source_object_content.endswith('.lnk'):
                destination_object_content_path = fr'{destination_object_path}\{Path(source_object_content).stem}'
                shortcut = shell.CreateShortCut(source_object_content_path)
                source_object_content_path = shortcut.Targetpath
            else:
                destination_object_content_path = fr'{destination_object_path}\{source_object_content}'
            if os.path.isdir(source_object_content_path):
                if not os.path.isdir(destination_object_content_path):
                    print('\nMaking directory ' + destination_object_content_path)
                    try: os.mkdir(destination_object_content_path)
                    except:
                        print('Make failed')
                        failed_actions.append(f'Make {destination_object_content_path}')
                recursive_selective_file_transfer(source_object_content_path, destination_object_content_path, objects_not_to_copy, paths_not_to_copy)
            else:
                if os.path.isfile(destination_object_content_path):
                    if replace_or_update_or_keep_or_forbid_existing_duplicate_files == 'r':
                        do_copy = 1
                    elif replace_or_update_or_keep_or_forbid_existing_duplicate_files == 'k':
                        do_copy = 0
                    elif replace_or_update_or_keep_or_forbid_existing_duplicate_files == 'u':
                        if os.stat(source_object_content_path).st_mtime > os.stat(destination_object_content_path).st_mtime:
                            do_copy = 1
                        else:
                            do_copy = 0
                    else:
                        raise Exception(f"The destination file, {destination_object_content_path}, already exists - to proceed, it must be deleted, renamed, or moved")
                else: do_copy = 1
                if do_copy == 1:
                    if not os.path.isdir(destination_object_path):
                        print('\nMaking directory ' + destination_object_path)
                        try: os.mkdir(destination_object_path)
                        except:
                            print('Make failed')
                            failed_actions.append(f'Make {destination_object_path}')
                    print(f'\nCopying {source_object_content_path} to {destination_object_content_path}')
                    try: shutil.copy2(source_object_content_path, destination_object_content_path)
                    except:
                        print('Transfer failed')
                        failed_actions.append(f'{source_object_content_path} to {destination_object_content_path}')
                        
    return failed_actions

shell = win32com.client.Dispatch("WScript.Shell")
for subdirectory_idx, source_and_destination_subdirectory in enumerate(source_and_destination_subdirectories):
    source_subdirectory = source_and_destination_subdirectory[0]
    destination_subdirectory = source_and_destination_subdirectory[1]
    source_subdirectory_path = fr'{top_source_directory_path}\{source_subdirectory}'
    destination_subdirectory_path = fr'{top_destination_directory_path}\{destination_subdirectory}'
    if os.path.exists(source_subdirectory_path) and source_subdirectory not in objects_not_to_copy:
        failed_actions = recursive_selective_file_transfer(source_subdirectory_path, destination_subdirectory_path, objects_not_to_copy, paths_not_to_copy)
        
et = time.time()
elapsed_time = et - st
if elapsed_time < 60:
    print('\nExecution time:', elapsed_time, 'seconds')
elif elapsed_time >= 60 and elapsed_time < 60*60:
    print('\nExecution time:', elapsed_time/60, 'minutes')
else:
    print('\nExecution time:', elapsed_time/60/60, 'hours')

if failed_actions:
    failed_actions_text = '\n'.join(map(str, failed_actions))
    print(f'\nFailed actions: \n{failed_actions_text}')