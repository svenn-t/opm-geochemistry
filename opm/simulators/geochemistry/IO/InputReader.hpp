/*
 * MIT License
 *
 * Copyright (C) 2025 Aksel Hiorth
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE.
*/
#ifndef INPUT_READER_DEF_H
#define INPUT_READER_DEF_H

#include <algorithm>
#include <fstream>
#include <istream>
#include <iostream>
#include <map>
#include <sstream>
#include <string>
#include <vector>

/**
* Helper class used to scan through an input stream, searching for pre-specified keywords.
*
* Currently, two types of keywords are supported:
*
* 1) Unique keywords that have *exactly* one associated value (a string):
*
* - If present, such a keyword is assumed to come first in a line in the input stream.
* - The corresponding value should follow immediately thereafter (separated by white space).
* - If a keyword is entered multiple times, the value at the last occurrence will be stored.
*
* 2) "Block" keywords that store more complex input as vectors of strings (for further processing):
*
* - Block keywords are defined by the combination of a "base" keyword and a unique name (OPTIONAL; default is the empty string).
* - Block keywords with the same "base keyword", but different names, are therefore treated as different.
* (concrete example: "solution 0" vs "solution 1")
* - Each block keyword is terminated by an appropriate closing tag (a string). Several such strings may be specified, in which
* case the reading of the keyword is terminated the first time one of them is encountered
* (example: "/end" or "/ end"; this way the code becomes more robust when it comes to user input errors)
* - As with the simple keywords: if a block keyword is entered multiple times, only the last occurrence will be stored.
*
* Other things to note:
*
* - Lines starting with either "#", "*", or "-" are assumed to be comment lines.
* - If case_sensitive_ == false (default behaviour), keywords have to be searched for in UPPER CASE.
* - Otherwise, keywords can be searched for with a mix of upper and lower case letters.
*/
class InputReader{

  public:

    explicit InputReader(const std::map<std::string, std::string>& default_key_value_pairs, bool case_sensitive=false);

    InputReader(const InputReader& prmReader) = delete;
    InputReader(const InputReader&& prmReader) = delete;
    InputReader& operator=(const InputReader& prmReader) = delete;

    void removeAsteriskAsCommentMarker();

    /**
     * Read (key,value) pairs from an input stream.
     */
    void read_key_value_pairs(std::istream& input_stream);

    /**
     * User-added simple keywords.
     */
    void define_simple_keyword(const std::string& keyword, const std::string& default_value);

    /**
     * @return An std::map with the (keyword, value)-pairs used.
     */
    const std::map<std::string, std::string>& get_key_value_pairs() const;

    /**
     * @return An std::map with all the recorded block keywords and their contents.
     */
    const std::map<std::string, std::vector<std::string>>& get_all_block_keywords_content() const;

    /**
     * @return Either the value at the searched-for key, or the empty string
     *         (if not present in the list of pre-defined keywords).
     */
    std::string get_simple_keyword_value(const std::string& key);

    /**
     * @return The contents for a specific block keyword (e.g., "Solution 2").
     *          If nothing is found, return an empty vector of strings.
     */
    std::vector<std::string>& get_block_keyword_content(const std::string& unique_block_key);

    /**
     * Defines the start and end of a block.
     */
    void define_block_keyword(const std::string& block_begin, const std::vector<std::string>& block_end);

    /**
     * Scans through the input stream and stores input for the pre-defined block
     * keywords (exclude commented lines).
     *
     * @return A number > 0 if not all pre-specified keywords could be read.
     */
    int read_block_keywords(std::istream& input_stream);

    /**
     * Defines the start and end of the main block used.
     * Useful when we are interested in scanning a subset of a larger file.
     */
    void define_main_block(const std::string& main_block_name, const std::vector<std::string>& main_block_end_keywords);

    /** 
     * A high level function to read all keywords, including block keywords.
     */
    void read(std::istream& input_stream);

    /**
     * @return The number of duplicate block keywords.
     */
    int size(const std::string& block_name);

    /**
     * @return The number of lines for the given unique block name.
     */
    int block_size(const std::string& unique_block_name);

    /**
     * @return The number of lines for ALL defined blocks.
     */
    int block_size();

    /**
     * @return A vector of unique names.
     */
    std::vector<std::string> unique_names(const std::string& block_name);

    /**
     * @return A vector of unique names (all the blocks read).
     */
    std::vector<std::string> get_block_names_read();

  private:

    bool verbose_output_;
    bool case_sensitive_;

    std::vector<std::string> comment_markers_;
    std::string value_separator_;// odn (3/2-2020): not used atm.

    // empty string reference to string is returned whenever block_content or simple_key not found
    std::vector<std::string> empty_block_content_;
    std::string empty_simple_key_value_;

    // Key-value pairs:
    std::map<std::string, std::string> default_params_;
    std::map<std::string, std::string> used_params_;

    // block keywords and their names, e.g.: "solution" --> {"solution 0", "solution 1", ..}
    std::map<std::string, std::vector<std::string>> blocks_names_;
    // block keywords and their allowed closing strings
    std::map<std::string, std::vector<std::string>> blocks_end_keywords_;
    // unique block identifiers (keyword+name) and associated content
    std::map<std::string, std::vector<std::string>> blocks_content_;

    // Main Block - if not set, the whole input stream is read, otherwise we only look between start and stop
    std::string main_block_name_;
    std::vector<std::string> main_block_end_keywords_;

    // Private methods
    bool is_comment(const std::string& line);
    void reset_key_value_pairs_to_default_values();
    bool in_main_block(std::istream& input_stream, std::string& line, int no_lines_read) const;
    void find_bof_eof(std::istream& input_stream);

    long bof_; // start of input stream
    long eof_; // end of input stream
    int  main_block_lines_; // number of lines to be read from bof_
    void rewind(std::istream& input_stream) const;
    bool end_of_main_block_defined_;

};

#endif
