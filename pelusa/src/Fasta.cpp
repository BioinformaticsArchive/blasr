#include "Fasta.h"
#include <iostream>


FastaRecord::FastaRecord()
{
}

bool FastaRecord::parseRecord(FILE* fasta_file)
{
	char name[MAX_LINE];    
 	char desc[MAX_LINE];    
  	char buffer[MAX_SEQUENCE];
  	
  	// Read the title line.
  	if (!readTitleLine(fasta_file, name, desc)) {
    	return false;
  	}
  
  	// Read the sequence.
  	if (!readRawSequence(fasta_file,  buffer)) {
    	fprintf(stderr, "Sequence %s is too long.\n", name);
    	exit(1);
  	}
	this->name = string(name);
	this->description = string(desc);
	this->sequence = string(buffer);
  	return true;
}


/**
 */
bool FastaRecord::readTitleLine
  (FILE* fasta_file,
   char* name,
   char* description)
{
    char id_line[MAX_LINE];  // Line containing the ID and comment.
    int a_char;                         // The most recently read character.

    // Read until the first occurrence of ">".
    while ((a_char = getc(fasta_file)) != '>') {
        // If we hit the end of the file, return false.
        if (a_char == EOF) {
            return false;
        }  
    }
    char* new_line = NULL;
    int line_length;
    size_t buf_length = 0;

    if((line_length =  getline(&new_line, &buf_length, fasta_file)) == -1){
        fprintf(stderr, "Error reading Fasta file.\n");
        exit(1);
    }
    strncpy(id_line, new_line, MAX_LINE-1);
    delete(new_line); // TODO this causes a mis matched free error

    // Remove EOL.
    id_line[strlen(id_line) - 1] = '\0';

    // Extract the ID from the beginning of the line.
    if (sscanf(id_line, "%s", name) != 1) {
        ; // die
    }

    // Store the rest of the line as the comment.
    strcpy(description, &(id_line[strlen(name)+1]));

    return true;
}


/****************************************************************************
 ****************************************************************************/
bool FastaRecord::readRawSequence
  (FILE* fasta_file,   // Input Fasta file.
   char* raw_sequence)
{
  int a_char;
  unsigned int i_seq = 0;
  bool return_value = true;

  // Read character by character.
  while ((a_char = getc(fasta_file)) != EOF) {

    // Check for the beginning of the next sequence.
    if (a_char == '>') {
      ungetc(a_char, fasta_file);
      break;
    }

    // Skip non-alphabetic characters.
    if ((a_char == ' ')     || 
        (a_char == '\t')    ||
        (a_char == '\n')    || 
        (a_char == '\r')) {
        ; // skip
    } else {
      a_char = toupper((int)a_char);
      raw_sequence[i_seq] = a_char;
      i_seq++;
    }
  }
  raw_sequence[i_seq] = '\0';

  return(return_value);
}
