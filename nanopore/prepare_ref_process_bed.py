import sys
import os

def process_bed(input_file, output_file):
    print(f"Processing input file: {input_file}")
    print(f"Output will be saved to: {output_file}")
    
    awk_command = f"""
    awk '
    BEGIN {{
        i = 0;
    }}
    {{  
        if (i % 2 == 0) {{
            # First line of the pair
            coordinate_1 = $2;
            name = $4; 
            sub(/_.*/, "", name);
        }} else {{
            # Second line of the pair
            coordinate_2 = $3;
            size_1 = coordinate_2 - coordinate_1;
            size_2 = size_1 - 250;
            size_3 = size_1 + 250;
            printf "%s\\t%d\\t%d\\t%d\\n", name, size_1, size_2, size_3 >> "{output_file}";
        }}
        i++;
    }}
    ' {input_file}
    """
    print(f"Running command: {awk_command}")
    os.system(awk_command)
    print("Command executed")

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python prepare_ref_process_bed.py <input_file> <output_file>")
        sys.exit(1)
    
    input_file = sys.argv[1]
    output_file = sys.argv[2]
    print(f"Arguments received: input_file={input_file}, output_file={output_file}")
    process_bed(input_file, output_file)
    print("Script finished")

