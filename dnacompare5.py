import tkinter as tk
from tkinter import ttk, filedialog, messagebox
from Bio import SeqIO
from Bio.SeqUtils import gc_fraction
from Bio import pairwise2

# Global variable to hold results for sorting and searching
results = []

def compute_similarity(seq1, seq2):
    """Calculate percentage identity between two DNA sequences using global alignment with penalties."""
    # Perform global alignment with penalties for mismatches and gaps
    alignments = pairwise2.align.globalms(seq1, seq2, match=1, mismatch=-1, open=-2, extend=-0.5)
    best_alignment = alignments[0]  # Take the best alignment
    matches = sum(1 for a, b in zip(best_alignment[0], best_alignment[1]) if a == b)
    aligned_length = len(best_alignment[0].replace('-', ''))  # Exclude gaps from alignment length
    return (matches / aligned_length) * 100

def classify_relationship(similarity):
    """Classify relationship based on similarity percentage."""
    if similarity > 90:
        return "Parent/Sibling"
    elif 65 <= similarity <= 90:
        return "Distant Relative"
    else:
        return "Stranger"

# Add a global variable to store the first person's details
first_person_id = ""
first_person_gc = 0.0

def load_fasta():
    """Load FASTA file and process DNA relationships."""
    global results, first_person_id, first_person_gc
    file_path = filedialog.askopenfilename(
        filetypes=[("FASTA files", "*.fasta *.fa"), ("All files", "*.*")]
    )
    if not file_path:
        return

    try:
        sequences = list(SeqIO.parse(file_path, "fasta"))
        if not sequences:
            messagebox.showerror("Error", "No sequences found in the FASTA file.")
            return

        # Extract details of the first sequence
        first_seq = str(sequences[0].seq)
        first_person_id = sequences[0].id
        first_person_gc = gc_fraction(sequences[0].seq) * 100

        results = []  # Clear previous results

        for i, record in enumerate(sequences[1:], start=1):
            similarity = compute_similarity(first_seq, str(record.seq))
            gc_content = gc_fraction(record.seq) * 100
            relationship = classify_relationship(similarity)
            results.append({
                "id": record.id,
                "similarity": similarity,
                "gc_content": gc_content,
                "relationship": relationship,
            })

        update_results()  # Display the results
    except Exception as e:
        messagebox.showerror("Error", f"An error occurred: {e}")

def update_results(sort_key=None, search_keyword=None):
    """Update the results display with optional sorting or filtering."""
    global first_person_id, first_person_gc
    filtered_results = results

    if search_keyword:
        search_keyword = search_keyword.lower()
        filtered_results = [
            result for result in results
            if search_keyword in result["id"].lower() or search_keyword in result["relationship"].lower()
        ]

    if sort_key:
        filtered_results = sorted(filtered_results, key=lambda x: x[sort_key], reverse=True)

    results_text.delete(1.0, tk.END)

    # Add the first person's details to the header
    results_text.insert(
        tk.END, 
        f"Analyzing relationships for {first_person_id} (GC Content: {first_person_gc:.2f}%):\n\n"
    )

    for i, result in enumerate(filtered_results, start=1):
        results_text.insert(
            tk.END,
            f"Sequence {i} ({result['id']}):\n"
            f"  Similarity: {result['similarity']:.2f}%\n"
            f"  GC Content: {result['gc_content']:.2f}%\n"
            f"  Relationship: {result['relationship']}\n\n",
        )

def clear_results():
    """Clear the results from the text box."""
    results_text.delete(1.0, tk.END)
    global results
    results = []

def on_sort_change(event):
    """Handle sort dropdown change."""
    sort_option = sort_var.get()
    sort_key_map = {
        "Similarity": "similarity",
        "GC Content": "gc_content",
        "Relationship": "relationship",
    }
    update_results(sort_key=sort_key_map.get(sort_option))

def search_results():
    """Filter results based on search keyword."""
    search_keyword = search_entry.get()
    if search_keyword:
        update_results(search_keyword=search_keyword)

# GUI Setup
root = tk.Tk()
root.title("DNA Relationship Analyzer (NOTE: RELATIONSHIP STATUS MAYBE BIASED, USE SIMILARITY PERCENTAGE FOR BETTER UNDERSTANDING)")
root.geometry("800x600")

frame = ttk.Frame(root, padding=10)
frame.pack(fill=tk.BOTH, expand=True)

# Top Buttons
button_frame = ttk.Frame(frame)
button_frame.pack(fill=tk.X, pady=5)

load_button = ttk.Button(button_frame, text="Load FASTA File", command=load_fasta)
load_button.pack(side=tk.LEFT, padx=5)

clear_button = ttk.Button(button_frame, text="Clear Results", command=clear_results)
clear_button.pack(side=tk.LEFT, padx=5)

# Sorting Dropdown
sort_var = tk.StringVar(value="Sort By")
sort_dropdown = ttk.OptionMenu(
    button_frame, sort_var, "Sort By", "Similarity", "GC Content", "Relationship", command=on_sort_change
)
sort_dropdown.pack(side=tk.RIGHT, padx=5)

# Search Bar
search_frame = ttk.Frame(frame)
search_frame.pack(fill=tk.X, pady=5)

search_label = ttk.Label(search_frame, text="Search:")
search_label.pack(side=tk.LEFT, padx=5)

search_entry = ttk.Entry(search_frame)
search_entry.pack(side=tk.LEFT, fill=tk.X, expand=True, padx=5)

search_button = ttk.Button(search_frame, text="Search", command=search_results)
search_button.pack(side=tk.LEFT, padx=5)

# Results Display
results_text = tk.Text(frame, wrap=tk.WORD, width=80, height=30, font=("Arial", 12))
results_text.pack(side=tk.LEFT, fill=tk.BOTH, expand=True, padx=5, pady=5)

scrollbar = ttk.Scrollbar(frame, command=results_text.yview)
results_text.config(yscrollcommand=scrollbar.set)
scrollbar.pack(side=tk.RIGHT, fill=tk.Y)

root.mainloop()
