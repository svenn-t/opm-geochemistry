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
#include <opm/simulators/geochemistry/IO/InputReader.hpp>

#include <opm/simulators/geochemistry/Utility/HelpfulStringMethods.hpp>

InputReader::InputReader(const std::map<std::string, std::string>& default_key_value_pairs, bool case_sensitive)
    : verbose_output_(false)
    , case_sensitive_(case_sensitive)
    , comment_markers_{ "#", "*", "-","\r","\n" }
    , value_separator_("\t")
    , empty_block_content_()
    , empty_simple_key_value_()
    , default_params_()
    , used_params_()
    , blocks_names_()
    , blocks_end_keywords_()
    , blocks_content_()
    , main_block_name_()
    , main_block_end_keywords_()
    , bof_(-1)
    , eof_(-1)
    , main_block_lines_(-1)
    , end_of_main_block_defined_(false)
{

    if(!case_sensitive_)
    {
        for (const auto& default_key_value_pair : default_key_value_pairs)
        {
            std::string key = default_key_value_pair.first;
            default_params_[to_upper_case(key)] = default_key_value_pair.second;
        }
    }
    else
    {
        default_params_ = default_key_value_pairs;
    }

    used_params_ = default_params_;

}

void  InputReader::removeAsteriskAsCommentMarker() {
    for (auto it = comment_markers_.begin(); it != comment_markers_.end(); ++it) {
        if ((*it) == "*") {
            comment_markers_.erase(it);
            return;
        }
    }
}

std::vector<std::string> InputReader::get_block_names_read()
{
    std::vector<std::string> names;
    for (const auto& block_pair : blocks_content_) names.push_back(block_pair.first);
    return names;
}

void InputReader::read_key_value_pairs(std::istream& input_stream)
{
    if (input_stream.fail())
    {
        std::cout << " InputReader: could not read input!\n";
        exit(1);
        return;
    }
    bof_ = input_stream.tellg(); // get position 
    reset_key_value_pairs_to_default_values();

    std::string line;
    int no_lines_read = 0;
    while (in_main_block(input_stream, line, no_lines_read))
    {
        no_lines_read++;
        if (line.empty() || is_comment(line))
        {
            continue;
        }
        else
        {
            std::vector<std::string> lineTab = split_string(line);

            if (lineTab.size() < 2)
            {
                continue;
            }

            // We (currently) assume a single value for each keyword:
            std::string key = lineTab[0];
            std::string value = lineTab[1];

            if(!case_sensitive_)
            {
                key = to_upper_case(key);
            }

            if (used_params_.find(key) != used_params_.end())
            {
                used_params_[key] = value;

            }
            else if(verbose_output_)
            {
                std::cout << "Keyword " << key << " was not found!\n";
            }
        }
    }
    rewind(input_stream);
}

const std::map<std::string, std::string>& InputReader::get_key_value_pairs() const
{
    return used_params_;
}

const std::map<std::string, std::vector<std::string>>& InputReader::get_all_block_keywords_content() const
{
    return blocks_content_;
}

void InputReader::define_main_block(const std::string& main_block_name, const std::vector<std::string>& main_block_end_keywords)
{
    main_block_name_ = case_sensitive_ ? main_block_name : to_upper_case(main_block_name); 
    
    if (!case_sensitive_)
    {
        main_block_end_keywords_.clear();
        for (auto it = main_block_end_keywords.cbegin(); it < main_block_end_keywords.cend(); ++it)
        {
            main_block_end_keywords_.push_back(to_upper_case(*it));
        }
    }
    else
        main_block_end_keywords_ = main_block_end_keywords;

    end_of_main_block_defined_ = true;
}

std::string InputReader::get_simple_keyword_value(const std::string& key)
{
    const std::string key_search = case_sensitive_ ? key : to_upper_case(key);
    if(used_params_.find(key_search) != used_params_.end())
    {
        return used_params_.at(key_search);
    }
    return empty_simple_key_value_;
}

bool InputReader::in_main_block(std::istream& input_stream, std::string& line, int no_lines_read) const
{
    if (input_stream.eof())
    {
        return false;
    }
    
    if (end_of_main_block_defined_)
    {
        if (no_lines_read < main_block_lines_ || main_block_lines_ == -1)
        {
            std::getline(input_stream, line);
            line = trim_left_right(line);
            return true;
        }
        else
        {
            return false;
        }
    }
    else
    {
        std::getline(input_stream, line);
    }
    return true;
}

std::vector<std::string>& InputReader::get_block_keyword_content(const std::string& unique_block_key)
{
    const std::string key_search = case_sensitive_ ? unique_block_key : to_upper_case(unique_block_key);      
    if(blocks_content_.find(key_search) != blocks_content_.end())
    {
        return blocks_content_.at(key_search);
    }
    return empty_block_content_;
}

void InputReader::define_block_keyword(const std::string& block_begin, const std::vector<std::string>& block_end)
{
    if(case_sensitive_)
    {
        blocks_names_[block_begin] = std::vector<std::string>();
        blocks_end_keywords_.insert(std::pair <std::string, std::vector<std::string>>(block_begin, block_end));
    }
    else
    {
        // Note: end keywords are still case-sensitive...
        std::vector<std::string> end_keywords;
        for (std::vector<std::string>::const_iterator it = block_end.begin(); it < block_end.end(); ++it)
        {
            end_keywords.push_back(to_upper_case(*it));
        }
        
        blocks_names_[to_upper_case(block_begin)] = std::vector<std::string>();
        blocks_end_keywords_.insert(std::pair <std::string,std::vector<std::string>>(to_upper_case(block_begin), end_keywords));
        
    }
}

void InputReader::define_simple_keyword(const std::string& keyword, const std::string& default_value)
{
    std::string KEYWORD;
    if (case_sensitive_) KEYWORD = keyword;
    else KEYWORD = to_upper_case(keyword);
    
    used_params_[KEYWORD] = default_value;
    default_params_[KEYWORD] = default_value;
}

int InputReader::read_block_keywords(std::istream& input_stream)
{
    // Continually updated as new blocks are read
    bool inside_block = false;
    bool reached_end_of_block = false;
    
    std::string current_block_key;
    std::string current_block_name;
    std::vector<std::string> current_block_content;
    std::string proposed_key;
    std::string current_line;
    
    bof_ = input_stream.tellg();
    int no_lines_read = 0;
    while (in_main_block(input_stream, current_line, no_lines_read))
    {
        ++no_lines_read;
        if (current_line.empty() || is_comment(current_line))
        {
            continue;
        }
        else
        {
            const std::vector<std::string> lineTab = split_string(current_line);
            proposed_key = "";
            if (!lineTab.empty())
                proposed_key = case_sensitive_ ? lineTab[0] : to_upper_case(lineTab[0]);
                        
            if (inside_block)
            {
                current_block_content.push_back(trim_left_right(current_line));
            }
            else if(blocks_end_keywords_.find(proposed_key) != blocks_end_keywords_.end())
            {
                // We have found a new 'block keyword'
                inside_block = true;
                current_block_key = proposed_key;
                current_block_name = lineTab.size() > 1 ? lineTab[1] : ""; // typically a number
                if (!case_sensitive_) current_block_name = to_upper_case(current_block_name);
            }

            if(inside_block)
            {
                // Check if the current line starts with one of the closing tags:
                const std::string line_to_match = case_sensitive_ ? current_line : to_upper_case(current_line);
                // ah 09.06.20 - better to check proposed key, more robust - 2021 12.11 change back.
                for(const auto& closing_tag: blocks_end_keywords_[current_block_key])
                {
                    if(!line_to_match.compare(0, closing_tag.size(), closing_tag))
                    {
                        reached_end_of_block = true;
                        break;
                    }
                }
                
                if(reached_end_of_block)
                {
                    // We are done reading the current keyword:
                    inside_block = false;
                    reached_end_of_block = false;

                    // Remove 'closing tag' of the block keyword, i.e., keep the actual content only:
                    current_block_content.erase(current_block_content.end() - 1);

                    // Store contents, and make sure to distinguish between different block names:
                    std::string unique_key = current_block_name.empty() ? current_block_key : (current_block_key + " " + current_block_name);
                    if(!case_sensitive_) unique_key = to_upper_case(unique_key);

                    blocks_names_.at(current_block_key).push_back(current_block_name);
                    blocks_content_[unique_key] = current_block_content;

                    // Reset in case we have to read more input:
                    current_block_key = "";
                    current_block_name = "";
                    current_block_content = std::vector<std::string>();
                }
            } // end inside block
        }
    }
    rewind(input_stream); // set back pointer before exiting
    if(inside_block && !reached_end_of_block)
    {
        std::cout << "Error when reading block keyword " << current_block_key << "...\n";
        return 1;
    }
    return 0;
}

void InputReader::read(std::istream& input_stream)
{
    if(input_stream.fail())
    {
        std::cout << "Error when trying to read input\n";
        exit(1);
    }
    find_bof_eof(input_stream);
    read_block_keywords(input_stream);
    read_key_value_pairs(input_stream);
//    input_stream.seekg(eof_); // bug??? ah: 16.06.2022
    input_stream.seekg(bof_);
}

int InputReader::size(const std::string& block_name)
{
    return blocks_names_[block_name].size();
}

int InputReader::block_size(const std::string& unique_block_name)
{
    return blocks_content_[unique_block_name].size();
}

int InputReader::block_size()
{
    std::size_t total_size = 0;
    for (const auto& block_pair : blocks_content_) total_size += block_pair.second.size();
    return total_size;
}

std::vector<std::string> InputReader::unique_names(const std::string& block_name)
{
    const std::string w = case_sensitive_ ? block_name : to_upper_case(block_name);
    
    std::vector<std::string> ww;
    for (auto it = blocks_names_[w].cbegin(); it < blocks_names_[w].cend(); ++it){
        ww.push_back(trim_left_right(block_name + " " + *it));
    }
    return ww;
}

// ************************************************** PRIVATE METHODS **************************************************
bool InputReader::is_comment(const std::string& line)
{
    for (auto it = comment_markers_.cbegin(); it < comment_markers_.cend(); ++it)
    {
        if (!line.empty() && line[0] == (*it)[0])
            return true;
    }
    return false;
}

void InputReader::reset_key_value_pairs_to_default_values()
{
    used_params_ = default_params_;
}

// Set back stream. If eof is reached, the stream must be cleared.
void InputReader::rewind(std::istream& input_stream) const
{
    if (input_stream.eof()) input_stream.clear();

    input_stream.seekg(bof_);
}

// Finds the beginning and end of the file stream.
void InputReader::find_bof_eof(std::istream& input_stream)
{
    input_stream.seekg(0, input_stream.end);
    eof_ = input_stream.tellg();
    input_stream.seekg(0, input_stream.beg);
    bof_ = input_stream.tellg();
    if (!end_of_main_block_defined_)
        return;
    else
    {
        std::string current_line, proposed_key, current_block_key;
        bool inside_block = false;
        bool reached_end_of_block = false;
        bool inside_main_block = false;
        bool reached_end_of_main_block = false;
        // search from beginning until the main block keyword
        while (!reached_end_of_main_block)
        {
            ++main_block_lines_;
            std::getline(input_stream, current_line);
            const std::vector<std::string> lineTab = split_string(current_line);
            proposed_key = "";
            if (!lineTab.empty())
                proposed_key = case_sensitive_ ? lineTab[0] : to_upper_case(lineTab[0]);
            
            if (proposed_key == main_block_name_)
            {
                inside_main_block = true;
                bof_ = input_stream.tellg();
                main_block_lines_ = 0; // reset counter and count from here
            }
            if (inside_main_block)
            {
                if (blocks_end_keywords_.find(proposed_key) != blocks_end_keywords_.end())
                {
                    // We found a new 'block keyword':
                    inside_block = true;
                    current_block_key = proposed_key;
                }

                if (inside_block)
                {           
                    // Check if the current line starts with one of the closing tags:
                    const std::string line_to_match = case_sensitive_ ? current_line : to_upper_case(current_line);
                    for (const auto& closing_tag : blocks_end_keywords_[current_block_key])
                    {
                        if (!line_to_match.compare(0, closing_tag.size(), closing_tag))
                        {
                            reached_end_of_block = true;
                            break;
                        }
                    }
                    if (reached_end_of_block)
                    {
                        // We are done reading the current keyword:
                        inside_block = false;
                        reached_end_of_block = false;

                        // Reset in case we have to read more input:
                        current_block_key = "";
                    }
                }
                else  // check if end of main block is reached
                {
                    const std::string line_to_match = case_sensitive_ ? current_line : to_upper_case(current_line);
                    for (const auto& closing_tag : main_block_end_keywords_)
                    {
                        if (!line_to_match.compare(0, closing_tag.size(), closing_tag))
                        {
                            reached_end_of_main_block = true;
                            eof_ = input_stream.tellg();
                            input_stream.clear();
                            input_stream.seekg(bof_);
                            break;
                        }
                    }
                }
            }
            if (input_stream.eof())
            {
                std::cout << " Error: main block " << main_block_name_ << " not correctly defined, maybe missing or wrong spelling of end-of-block\n";
                input_stream.seekg(bof_);
                return;
            }
        }
    }
}

