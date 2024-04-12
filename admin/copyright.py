#!/usr/bin/env python3
#
# Membrane Dynamics in 3D using Discrete Differential Geometry (Mem3DG).
#
# Copyright 2024- The Mem3DG Authors
# and the project initiators Cuncheng Zhu, Christopher T. Lee, and
# Padmini Rangamani.
#
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at http://mozilla.org/MPL/2.0/.
#
# Please help us support Mem3DG development by citing the research
# papers on the package. Check out https://github.com/RangamaniLabUCSD/Mem3DG/
# for more information.

"""Checks and/or updates copyright headers.

This script is shamelessly borrowed from the GROMACS project and adapted for use with Mem3DG.
"""

import datetime
import os.path
import re
import sys

from optparse import OptionParser


class CopyrightState(object):
    """Information about an existing (or non-existing) copyright header."""

    def __init__(self, has_copyright, is_correct, is_newstyle, first_year):
        self.has_copyright = has_copyright
        self.is_correct = is_correct
        self.is_newstyle = is_newstyle
        self.first_year = first_year


class CopyrightChecker(object):
    """Logic for analyzing existing copyright headers and generating new ones."""

    _header = [
        "",
        "Membrane Dynamics in 3D using Discrete Differential Geometry (Mem3DG).",
        "",
    ]
    _copyright = "Copyright {0}- The Mem3DG Authors"
    _footer = """
and the project initiators Cuncheng Zhu, Christopher T. Lee, and
Padmini Rangamani.

This Source Code Form is subject to the terms of the Mozilla Public
License, v. 2.0. If a copy of the MPL was not distributed with this
file, You can obtain one at http://mozilla.org/MPL/2.0/.

Please help us support Mem3DG development by citing the research
papers on the package. Check out https://github.com/RangamaniLabUCSD/Mem3DG/
for more information.
""".strip().splitlines()

    def check_copyright(self, comment_block):
        """Analyze existing copyright header for correctness and extract information."""
        copyright_re_new = r"Copyright ([0-9]{4})- The Mem3DG Authors"
        copyright_re_old = r"Copyright \(c\) (([0-9]{4})):?.*"
        has_copyright = False
        is_newstyle = True
        is_correct = True
        found_old_style_line = False
        next_header_line = 0
        next_footer_line = 0
        first_year = ""
        for line in comment_block:
            if "Copyright" in line:
                has_copyright = True
                old_match = re.match(copyright_re_old, line)
                match = re.match(copyright_re_new, line)
                if (old_match and not match) and not found_old_style_line:
                    is_newstyle = False
                    found_old_style_line = True
                    get_first_year_re = r"Copyright \(c\) ([0-9]{4}[,-:]).*"
                    first_year_line_match = re.match(get_first_year_re, line)
                    first_year_line = first_year_line_match.group(1)
                    # decide which separator to use for finding first year
                    if ":" in first_year_line:
                        pos = first_year_line.find(":")
                        first_year = first_year_line[:pos]
                    else:
                        raise Exception("Can't find expected separator!")
                elif match:
                    first_year = match.group(1)
                    new_line = self._copyright.format(first_year)
                    if line != new_line:
                        is_correct = False
                if next_header_line != -1 or next_footer_line != 0:
                    is_correct = False
                continue
            if next_header_line >= 0:
                if line == self._header[next_header_line]:
                    next_header_line += 1
                    if next_header_line >= len(self._header):
                        next_header_line = -1
                else:
                    is_correct = False
                    is_newstyle = False
            elif next_footer_line >= 0:
                if line == self._footer[next_footer_line]:
                    next_footer_line += 1
                    if next_footer_line >= len(self._footer):
                        next_footer_line = -1
                else:
                    is_correct = False
            else:
                is_correct = False
        if next_header_line != -1 or next_footer_line != -1:
            is_correct = False

        return CopyrightState(has_copyright, is_correct, is_newstyle, first_year)

    def process_copyright(self, state, options, first_year, current_year, reporter):
        """Determine whether a copyright header needs to be updated and report issues."""
        need_update = False

        if state.first_year:
            first_year = state.first_year
            if not current_year >= first_year:
                if options.update_year:
                    need_update = True
                    first_year = current_year
                if options.check or not need_update:
                    reporter.report("initial copyright year wrong")
                else:
                    reporter.report("initial copyright year corrected")

        if not state.has_copyright:
            if options.add_missing:
                need_update = True
                first_year = current_year
            if options.check or not need_update:
                reporter.report("copyright header missing")
            elif options.add_missing:
                reporter.report("copyright header added")
        else:
            if not state.is_newstyle:
                if options.replace_header:
                    need_update = True
                if options.check or not need_update:
                    reporter.report("copyright header incorrect")
                else:
                    reporter.report("copyright header replaced")
            elif not state.is_correct:
                if options.update_header:
                    need_update = True
                if options.check or not need_update:
                    reporter.report("copyright header outdated")
                else:
                    reporter.report("copyright header updated")

        return need_update, first_year

    def get_copyright_text(self, first_year):
        """Construct a new copyright header."""
        output = []
        output.extend(self._header)

        output.append(self._copyright.format(first_year))
        output.extend(self._footer)
        return output


class Reporter(object):
    """Wrapper for reporting issues in a file."""

    def __init__(self, reportfile, filename):
        self._reportfile = reportfile
        self._filename = filename

    def report(self, text):
        self._reportfile.write(self._filename + ": " + text + "\n")


class CommentHandlerC(object):
    """Handler for extracting and creating C-style comments."""

    def extract_first_comment_block(self, content_lines):
        if not content_lines or not (
            content_lines[0].startswith("/*") or content_lines[0].startswith("//")
        ):
            return ([], 0)
        comment_block = [content_lines[0][2:].strip()]
        line_index = 1
        while line_index < len(content_lines):
            line = content_lines[line_index]
            # Detect end of block
            if (
                "*/" in content_lines[line_index]
                or len(content_lines[line_index]) == 0
                or "#include" in content_lines[line_index]
            ):
                break
            comment_block.append(line.lstrip("* ").lstrip("// ").rstrip())
            line_index += 1
        return (comment_block, line_index + 1)

    def create_comment_block(self, lines):
        output = []
        output.append(("/* " + lines[0]).rstrip())
        output.extend([(" * " + x).rstrip() for x in lines[1:]])
        output.append(" */")
        return output


class CommentHandlerSimple(object):
    """Handler for extracting and creating sh-style comments.

    Also other comments of the same type, but with a different comment
    character are supported."""

    def __init__(self, comment_char):
        self._comment_char = comment_char

    def extract_first_comment_block(self, content_lines):
        if not content_lines or not content_lines[0].startswith(self._comment_char):
            return ([], 0)
        comment_block = []
        line_index = 0
        while line_index < len(content_lines):
            line = content_lines[line_index]
            if not line.startswith(self._comment_char):
                break
            comment_block.append(line.lstrip(self._comment_char + " ").rstrip())
            line_index += 1
        while line_index < len(content_lines):
            line = content_lines[line_index].rstrip()
            if len(line) > 0 and line != self._comment_char:
                break
            line_index += 1
        return (comment_block, line_index)

    def create_comment_block(self, lines):
        output = []
        output.extend([(self._comment_char + " " + x).rstrip() for x in lines])
        output.append("")
        return output


comment_handlers = {
    "c": CommentHandlerC(),
    "tex": CommentHandlerSimple("%"),
    "sh": CommentHandlerSimple("#"),
}


def select_comment_handler(override, filename):
    """Select comment handler for a file based on file name and input options."""
    filetype = override
    if not filetype and filename != "-":
        basename = os.path.basename(filename)
        root, ext = os.path.splitext(basename)
        if ext == ".cmakein":
            dummy, ext2 = os.path.splitext(root)
            if ext2:
                ext = ext2
        if ext in (
            ".c",
            ".cpp",
            ".h",
            ".hpp",
        ):
            filetype = "c"
        elif ext in (".tex",):
            filetype = "tex"
        elif basename in ("CMakeLists.txt") or ext in (
            ".cmake",
            ".cmakein",
            ".in",
            ".py",
            ".sh",
            ".bash",
            ".csh",
            ".zsh",
        ):
            filetype = "sh"
    if filetype in comment_handlers:
        return comment_handlers[filetype]
    if filetype:
        sys.stderr.write("Unsupported input format: {0}\n".format(filetype))
    elif filename != "-":
        sys.stderr.write("Unsupported input format: {0}\n".format(filename))
    else:
        sys.stderr.write("No file name or file type provided.\n")
    sys.exit(1)


def create_copyright_header(current_year, language="c"):
    if language not in comment_handlers:
        sys.stderr.write("Unsupported language: {0}\n".format(language))
        sys.exit(1)
    copyright_checker = CopyrightChecker()
    comment_handler = comment_handlers[language]
    copyright_lines = copyright_checker.get_copyright_text(current_year)
    comment_lines = comment_handler.create_comment_block(copyright_lines)
    return "\n".join(comment_lines) + "\n"


def process_options():
    """Process input options."""
    parser = OptionParser()
    parser.add_option("-l", "--lang", help="Comment type to use (c or sh)")
    parser.add_option(
        "--first_year", help="Value for initial year to set in copyright file."
    )
    parser.add_option("-F", "--files", help="File to read list of files from")
    parser.add_option(
        "--check",
        action="store_true",
        help="Do not modify the files, only check the copyright (default action). "
        + "If specified together with --update, do the modifications "
        + "but produce output as if only --check was provided.",
    )
    parser.add_option(
        "--update-year",
        action="store_true",
        help="Update the copyright year if outdated",
    )
    parser.add_option(
        "--update-header",
        action="store_true",
        help="Update the copyright header if outdated",
    )
    parser.add_option(
        "--replace-header",
        action="store_true",
        help="Replace any copyright header with the current one",
    )
    parser.add_option(
        "--remove-old-copyrights",
        action="store_true",
        help="Remove copyright statements not in the new format",
    )
    parser.add_option(
        "--add-missing", action="store_true", help="Add missing copyright headers"
    )
    options, args = parser.parse_args()

    filenames = args
    if options.files:
        with open(options.files, "r") as filelist:
            filenames = [x.strip() for x in filelist.read().splitlines()]
    elif not filenames:
        filenames = ["-"]

    # Default is --check if nothing provided.
    if (
        not options.check
        and not options.update_year
        and not options.update_header
        and not options.replace_header
        and not options.add_missing
    ):
        options.check = True

    return options, filenames


def main():
    """Do processing as a stand-alone script."""
    options, filenames = process_options()
    first_year = options.first_year
    current_year = str(datetime.date.today().year)

    checker = CopyrightChecker()

    # Process each input file in turn.
    for filename in filenames:
        comment_handler = select_comment_handler(options.lang, filename)

        # Read the input file.  We are doing an in-place operation, so can't
        # operate in pass-through mode.
        if filename == "-":
            contents = sys.stdin.read().splitlines()
            reporter = Reporter(sys.stderr, "<stdin>")
        else:
            with open(filename, "r", encoding="utf-8") as inputfile:
                contents = inputfile.read().splitlines()
            reporter = Reporter(sys.stdout, filename)

        output = []
        # Keep lines that must be at the beginning of the file and skip them in
        # the check.
        if contents and (
            contents[0].startswith("#!/")
            or contents[0].startswith("%code requires")
            or contents[0].startswith("/* #if")
        ):
            output.append(contents[0])
            contents = contents[1:]
        # Remove and skip empty lines at the beginning.
        while contents and len(contents[0]) == 0:
            contents = contents[1:]

        # Analyze the first comment block in the file.
        comment_block, line_count = comment_handler.extract_first_comment_block(
            contents
        )

        state = checker.check_copyright(comment_block)

        need_update, first_year = checker.process_copyright(
            state, options, first_year, current_year, reporter
        )

        if options.remove_old_copyrights:
            need_update = True
            reporter.report("old copyrights removed")

        if need_update:
            # Remove the original comment if it was a copyright comment.
            if state.has_copyright:
                contents = contents[line_count:]
            new_block = checker.get_copyright_text(first_year)
            output.extend(comment_handler.create_comment_block(new_block))

        # Write the output file if required.
        if need_update or filename == "-":
            # Append the rest of the input file as it was.
            output.extend(contents)
            output = "\n".join(output) + "\n"
            if filename == "-":
                sys.stdout.write(output)
            else:
                with open(filename, "w", encoding="utf-8") as outputfile:
                    outputfile.write(output)


if __name__ == "__main__":
    main()
