#!/usr/bin/env python3
"""Developer-only registration inventory; names/arity are NOT functional parity.

Read the pinned adjacent C++ tree, model its Config registration constructor,
and compare name coverage with Rust REGISTERED_OPTIONS. No builds are run.

    python3 scripts/audit_command_options.py --write translation/command_options.json
    python3 scripts/audit_command_options.py --check translation/command_options.json
    python3 scripts/audit_command_options.py --probe .tmp/registry-probes.json

Only EXTRA and WITH_DNA branches occur in the registration span. Unexpected
syntax, unknown storage types, or new preprocessor conditions fail closed.
Probes always end with an invalid sentinel option, before Config/store returns;
CPU, address-space, output-size, core and wall-clock limits bound the children.
"""
from __future__ import annotations

import argparse
import ast
from collections import Counter
import hashlib
import json
from pathlib import Path
import re
import subprocess
import sys

VARIANTS = {
    "baseline": set(),
    "with_dna": {"WITH_DNA"},
    "extra": {"EXTRA"},
    "extra_with_dna": {"EXTRA", "WITH_DNA"},
}
TOKEN = re.compile(r'"(?:\\.|[^"\\])*"|\'(?:\\.|[^\'\\])*\'|[A-Za-z_]\w*|\d+(?:\.\d*)?(?:[eE][+-]?\d+)?[A-Za-z_]*|::|[^\s]')
COMMENTS = re.compile(r'"(?:\\.|[^"\\])*"|\'(?:\\.|[^\'\\])*\'|//[^\n]*|/\*.*?\*/', re.S)
TYPE_PREFIX = re.compile(r'^\s*((?:Option\s*<\s*(?:std::)?vector\s*<\s*(?:std::)?string\s*>\s*>|Option\s*<\s*\w+\s*>|(?:std::)?vector\s*<\s*(?:std::)?string\s*>|std::string|bool|double|float|int64_t|uint64_t|uint32_t|uint16_t|uint8_t|size_t|unsigned|int|long|Loc|string_vector|string))\s+(.+?)\s*;\s*$')


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def uncomment(text):
    def replace(match):
        value = match.group()
        return ''.join('\n' if c == '\n' else ' ' for c in value) if value.startswith(('//', '/*')) else value
    return COMMENTS.sub(replace, text)


def preprocess(text, definitions):
    """Blank inactive lines, preserving source offsets and line numbers."""
    output, stack, active = [], [], True
    for line in text.splitlines(keepends=True):
        directive = re.match(r'^\s*#\s*(\w+)\s*(.*?)\s*$', line)
        if directive:
            op, value = directive.groups()
            if op in {'ifdef', 'ifndef'}:
                if value not in {'EXTRA', 'WITH_DNA'}:
                    raise ValueError(f'Unknown registration preprocessor condition: {line.strip()}')
                condition = (value in definitions) == (op == 'ifdef')
                stack.append((active, condition, False))
                active = active and condition
            elif op == 'else':
                if not stack or stack[-1][2]:
                    raise ValueError('Unbalanced/repeated preprocessor else')
                parent, condition, _ = stack[-1]
                stack[-1] = (parent, condition, True)
                active = parent and not condition
            elif op == 'endif':
                if not stack:
                    raise ValueError('Unbalanced preprocessor endif')
                active = stack.pop()[0]
            else:
                raise ValueError(f'Unexpected registration directive: {line.strip()}')
            output.append(''.join('\n' if c == '\n' else ' ' for c in line))
        elif active:
            output.append(line)
        else:
            output.append(''.join('\n' if c == '\n' else ' ' for c in line))
    if stack:
        raise ValueError('Unterminated registration preprocessor block')
    return ''.join(output)


def arguments(text, opening):
    """Split one C++ call/list, keeping nested template/default expressions."""
    tokens = list(TOKEN.finditer(text, opening))
    if not tokens or tokens[0].start() != opening or tokens[0].group() not in '({':
        raise ValueError('Expected call or brace list')
    closing = ')' if tokens[0].group() == '(' else '}'
    stack = [closing]
    start, result = tokens[0].end(), []
    pairs = {'(': ')', '{': '}', '[': ']', '<': '>'}
    for token in tokens[1:]:
        value = token.group()
        if value in pairs:
            stack.append(pairs[value])
        elif value in {')', '}', ']', '>'}:
            if not stack or stack[-1] != value:
                raise ValueError(f'Unbalanced argument delimiter: {value}')
            stack.pop()
            if not stack:
                tail = text[start:token.start()].strip()
                if tail or result:
                    result.append(tail)
                return result, token.end()
        elif value == ',' and len(stack) == 1:
            result.append(text[start:token.start()].strip())
            start = token.end()
    raise ValueError('Unterminated argument list')


def string_literal(expression):
    strings = re.findall(r'"(?:\\.|[^"\\])*"', expression)
    if not strings or re.sub(r'"(?:\\.|[^"\\])*"', '', expression).strip():
        raise ValueError(f'Unsupported string expression: {expression}')
    return ''.join(ast.literal_eval(value) for value in strings)


def declarations(text, source, line_offset=0):
    found = {}
    for i, line in enumerate(uncomment(text).splitlines(), 1):
        match = TYPE_PREFIX.match(line)
        if not match:
            continue
        typename, names = match.groups()
        typename = re.sub(r'\s+', '', typename)
        for name in names.split(','):
            name = name.strip()
            if not re.fullmatch(r'\w+', name):
                continue  # Never infer declarators, initializers, or methods.
            found[name] = {'type': typename, 'source': source, 'line': i + line_offset}
    return found


def type_arity(declaration, min_expression):
    typename = declaration['type']
    if typename == 'bool':
        return {'kind': 'bool', 'minimum': 0, 'maximum': 0}
    if typename in {'string_vector', 'std::vector<string>', 'vector<string>', 'Option<string_vector>', 'Option<std::vector<std::string>>', 'Option<vector<string>>'}:
        if not re.fullmatch(r'\d+', min_expression):
            raise ValueError(f'Unknown vector minimum expression: {min_expression}')
        return {'kind': 'vector', 'minimum': int(min_expression), 'maximum': None}
    return {'kind': 'scalar', 'minimum': 1, 'maximum': 1}


def extract(repo):
    parser_header = (repo / 'diamond/src/util/command_line_parser.h').read_text(errors='replace')
    normalized = re.sub(r'\s+', '', uncomment(parser_header))
    parameter = '(conststd::vector<std::string>&v,constintmin_count)'
    signatures = {
        'scalar': 'template<typenameT>inlineboolcheck_pcount' + parameter + '{returnv.size()==1;}',
        'bool': 'template<>inlineboolcheck_pcount<bool>' + parameter + '{returnv.size()==0;}',
        'vector': 'template<>inlineboolcheck_pcount<std::vector<std::string>>' + parameter + '{returnv.size()>=(size_t)min_count;}',
        'optional_vector': 'template<>inlineboolcheck_pcount<Option<std::vector<std::string>>>' + parameter + '{returnv.size()>=(size_t)min_count;}',
    }
    for kind, signature in signatures.items():
        if signature not in normalized:
            raise ValueError(f'Unknown source parameter-count rule: {kind}')
    if 'map_short_[short_id]=o;' not in normalized:
        raise ValueError('Unknown short-option registration/shadow rule')
    parser_source = (repo / 'diamond/src/util/command_line_parser.cpp').read_text(errors='replace')
    documentation_source = parser_source.split('void CommandLineParser::print_documentation', 1)[1].split('void CommandLineParser::require', 1)[0]
    if not re.search(r'const size_t col1_width\s*=\s*25\s*;', documentation_source):
        raise ValueError('Unknown documentation column width')
    config_path = repo / 'diamond/src/basic/config.cpp'
    raw = config_path.read_text()
    begin = re.search(r'Config::Config\(int argc, const char \*\*argv, bool check_io, CommandLineParser& parser\)\s*\{', raw)
    if not begin:
        raise ValueError('Config registration constructor signature not found')
    end = raw.index('parser.store(argc, argv, command);', begin.end())
    registration = raw[begin.end():end]
    prefix_lines = raw[:begin.end()].count('\n')
    raw_header = (repo / 'diamond/src/basic/config.h').read_text()
    commands, groups, options, builds = {}, {}, {}, {}
    failures = []
    for variant, definitions in VARIANTS.items():
        code = preprocess(uncomment(registration), definitions)
        # All storage fields precede Config() and use only WITH_DNA.
        # Include fields after the command enum, without parsing methods.
        header = raw_header[raw_header.index('struct Config\n'):raw_header.index('\n\tConfig() {}')]
        header_lines = raw_header[:raw_header.index('struct Config\n')].count('\n')
        types = declarations(preprocess(uncomment(header), definitions), 'diamond/src/basic/config.h', header_lines)
        types.update(declarations(code, 'diamond/src/basic/config.cpp', prefix_lines))
        command_names, option_names, group_defs, short_map, short_history = [], [], {}, {}, {}
        line_at = lambda p: prefix_lines + code[:p].count('\n') + 1
        for match in re.finditer(r'\.add_command\s*\(', code):
            args, _ = arguments(code, code.index('(', match.start()))
            if len(args) != 3:
                raise ValueError('Unknown add_command signature')
            name, description, enum = string_literal(args[0]), string_literal(args[1]), args[2]
            value = {'name': name, 'description': description, 'enum_expression': enum, 'source_line': line_at(match.start()), 'variants': []}
            existing = commands.setdefault(name, value)
            if {k: v for k, v in existing.items() if k != 'variants'} != {k: v for k, v in value.items() if k != 'variants'}:
                raise ValueError(f'Inconsistent command definition: {name}')
            existing['variants'].append(variant)
            command_names.append(name)
        enum_to_name = {commands[name]['enum_expression']: name for name in command_names}
        for match in re.finditer(r'auto\s*&\s*(\w+)\s*=\s*parser\.add_group\s*\(', code):
            variable = match.group(1)
            args, _ = arguments(code, code.index('(', match.start()))
            if len(args) not in {2, 3}:
                raise ValueError('Unknown add_group signature')
            members, _ = arguments(args[1], 0)
            disabled = args[2] if len(args) == 3 else 'false'
            if disabled not in {'true', 'false'}:
                raise ValueError(f'Unknown group disabled expression: {disabled}')
            group = {'variable': variable, 'title': string_literal(args[0]), 'source_line': line_at(match.start()), 'command_enum_expressions': members, 'disabled_expression': disabled}
            group_defs[variable] = group
            common = groups.setdefault(variable, {'variable': variable, 'variants': {}})
            common['variants'][variant] = group
        for match in re.finditer(r'\b(\w+)\.add\s*\(\s*\)', code):
            variable = match.group(1)
            if variable not in group_defs:
                raise ValueError(f'Unregistered options group: {variable}')
            group = group_defs[variable]
            position = match.end()
            while True:
                while position < len(code) and code[position].isspace():
                    position += 1
                if code[position:position + 1] == ';':
                    break
                if code[position:position + 1] != '(':
                    raise ValueError(f'Unexpected option chain near line {line_at(position)}')
                opening = position
                args, position = arguments(code, position)
                if not 4 <= len(args) <= 6:
                    raise ValueError(f'Unknown option registration signature at {line_at(opening)}')
                name = string_literal(args[0])
                short = None if args[1] == '0' else ast.literal_eval(args[1])
                if short is not None and (not isinstance(short, str) or len(short) != 1):
                    raise ValueError(f'Invalid short option: {args[1]}')
                storage = args[3]
                declaration = types.get(storage)
                minimum = args[5] if len(args) > 5 else '1'
                if declaration is None:
                    arity = {'kind': 'unknown', 'reason': f'No recognized storage declaration: {storage}'}
                    failures.append(f'{variant}: --{name}: {arity["reason"]}')
                else:
                    arity = type_arity(declaration, minimum)
                value = {'name': name, 'short_id': short, 'short_expression': args[1], 'description': None if args[2] == 'format_str.str().c_str()' else string_literal(args[2]), 'description_expression': args[2], 'storage_expression': storage, 'storage_declaration': declaration, 'default_expression': args[4] if len(args) > 4 else 'T()', 'default_is_explicit': len(args) > 4, 'min_count_expression': minimum, 'min_count_is_explicit': len(args) > 5, 'arity': arity, 'group': variable, 'source_line': line_at(opening), 'variants': []}
                existing = options.setdefault(name, value)
                if {k: v for k, v in existing.items() if k != 'variants'} != {k: v for k, v in value.items() if k != 'variants'}:
                    raise ValueError(f'Inconsistent option registration: {name}')
                existing['variants'].append(variant)
                option_names.append(name)
                if short:
                    short_history.setdefault(short, []).append(name)
                    short_map[short] = name
        if len(option_names) != len(set(option_names)):
            raise ValueError(f'Duplicate long option IDs in {variant}: require explicit shadow modeling')
        builds[variant] = {'defines': sorted(definitions), 'workflow_permission_checks': 'EXTRA' not in definitions, 'command_count': len(command_names), 'commands': command_names, 'option_count': len(option_names), 'options': option_names, 'enabled_option_count': sum(group_defs[options[name]['group']]['disabled_expression'] == 'false' for name in option_names), 'disabled_options': [name for name in option_names if group_defs[options[name]['group']]['disabled_expression'] == 'true'], 'short_option_map': short_map, 'null_short_slot_final_option': option_names[-1], 'short_shadow_chains': {k: v for k, v in short_history.items() if len(v) > 1}, 'unregistered_group_command_enums': sorted({enum for g in group_defs.values() for enum in g['command_enum_expressions'] if enum not in enum_to_name})}
        builds[variant]['help_commands'] = [name for name in command_names if commands[name]['description']]
        builds[variant]['documentation_options'] = {
            name: [option for option in option_names
                   if commands[name]['enum_expression'] in group_defs[options[option]['group']]['command_enum_expressions']
                   and (options[option]['description'] or options[option]['description_expression'] == 'format_str.str().c_str()')]
            for name in command_names
        }
    rust = (repo / 'src/basic/config.rs').read_text()
    constant = re.search(r'pub const REGISTERED_OPTIONS: &str = ("(?:\\.|[^"\\])*");', rust)
    if constant is None:
        raise ValueError('Rust REGISTERED_OPTIONS constant not found')
    rust_names = ast.literal_eval(constant.group(1)).split()
    source_names = set(options)
    baseline = builds['baseline']
    if baseline['command_count'] != 21 or builds['extra']['command_count'] != 45:
        failures.append('Expected pinned inventory of 21 baseline and 24 additional EXTRA commands changed')
    if baseline['short_option_map'].get('w') != 'seed-weight':
        failures.append('Expected disabled seed-weight to shadow -w window changed')
    coverage = {'source_unique_names': len(source_names), 'rust_names': len(rust_names), 'missing_in_rust': sorted(source_names - set(rust_names)), 'rust_only_names': sorted(set(rust_names) - source_names), 'duplicate_rust_names': sorted(k for k, n in Counter(rust_names).items() if n > 1), 'meaning': 'Name coverage only; no claim of implemented command semantics, dispatch, arity, or output parity.'}
    if coverage['missing_in_rust'] or coverage['rust_only_names'] or coverage['duplicate_rust_names']:
        failures.append('Rust option-name coverage differs from source registration')
    artifact = {'schema_version': 1, 'source_commit': '1d162b4fefb5e6b4d868c24e6ad551cad2dcf246', 'purpose': 'Developer-only command/option registration inventory. This is not functionality or translation-completion proof.', 'source_sha256': {name: sha(repo / name) for name in ['diamond/src/basic/config.cpp', 'diamond/src/basic/config.h', 'diamond/src/util/command_line_parser.cpp', 'diamond/src/util/command_line_parser.h', 'src/basic/config.rs']}, 'registration_scope': {'constructor_begin_line': raw[:begin.start()].count('\n') + 1, 'parser_store_line': raw[:end].count('\n') + 1}, 'build_variants': builds, 'commands': list(commands.values()), 'groups': list(groups.values()), 'options': list(options.values()), 'rust_name_coverage': coverage, 'audit_failures': sorted(set(failures)), 'limitations': ['This targeted parser is not a general C++ compiler or full preprocessor; it fails on unknown registration syntax and conditions.', 'Types and arity are inferred only from recognized field/local declarations and CommandLineParser template rules. Post-store Config normalization, filesystem checks and command semantics are outside scope.', 'Registration defaults are preserved as expressions; implicit T() for Option<T> does not mean assignment because set_option_default(Option<T>) is a no-op.', 'The null short-ID slot is overwritten too but cannot name a regular argv option; the usable short-option map omits that null key.', 'Empty descriptions hide options/commands from help independently of registration; print_documentation itself does not test group.disabled.', 'All hidden options enter both long and short maps even when disabled. Later short IDs overwrite earlier IDs, including disabled seed-weight replacing -w window.', 'EXTRA disables workflow-membership checks; empty hidden group membership therefore permits enabled EXTRA hidden options on all registered commands.', 'WITH_DNA options register independently of EXTRA. The blastn command itself is registered only under EXTRA.', 'The dynamic outfmt description expression is retained without rendering its field table; known nonempty stream content establishes documentation visibility.', 'The source baseline binary validates only its compiled feature variant; alternate variants are statically modeled, not executed.', 'Source documentation pads option IDs to width 25 without a separator for longer IDs; a whitespace split cannot recover those names.', 'Bounded oracle probing requires POSIX Python resource limits; alternate build variants are not compiled by this tool.', 'CLI parser acceptance is not evidence that a command/option is implemented natively in Rust.']}
    return artifact



def probe_registry(repo, artifact, binary, destination):
    """Exercise parsing only, with an abort sentinel before Config returns."""
    import resource
    import tempfile

    sentinel = '__diamond_registry_audit_stop_79b1a__'
    if any(o['name'] == sentinel for o in artifact['options']):
        raise ValueError('Audit sentinel unexpectedly registered')
    binary = binary.resolve()
    expected_hash = 'd8d08608b96706db62e44036fc979f9a658e72f59fc19d8727546daa616c42d8'
    before = sha(binary)
    if before != expected_hash:
        raise ValueError('Oracle hash differs from immutable pinned baseline; do not infer its build flags')
    baseline = artifact['build_variants']['baseline']
    options = {o['name']: o for o in artifact['options']}
    groups = {g['variable']: g['variants']['baseline'] for g in artifact['groups']}
    commands = {c['name']: c for c in artifact['commands']}
    cases = []
    cap = 65536

    def limits():
        resource.setrlimit(resource.RLIMIT_CORE, (0, 0))
        resource.setrlimit(resource.RLIMIT_CPU, (2, 2))
        resource.setrlimit(resource.RLIMIT_AS, (512 << 20, 512 << 20))
        resource.setrlimit(resource.RLIMIT_FSIZE, (cap, cap))

    def run(kind, argv, expected_suffix, documentation=None):
        full = [str(binary), *argv]
        with tempfile.TemporaryFile() as stdout, tempfile.TemporaryFile() as stderr:
            try:
                result = subprocess.run(full, cwd=repo, stdout=stdout, stderr=stderr,
                                        preexec_fn=limits, timeout=3, check=False)
                stdout.seek(0)
                stderr.seek(0)
                out = stdout.read(cap + 1).decode(errors='replace')
                err = stderr.read(cap + 1).decode(errors='replace')
                passed = result.returncode == 1 and err.endswith(expected_suffix + '\n')
                record = {'kind': kind, 'args': argv, 'status': result.returncode,
                          'expected_error': expected_suffix, 'stdout': out, 'stderr': err}
                if documentation is not None:
                    option_lines = [line for line in out.splitlines() if line.startswith('--')]
                    actual_names, rendering_errors = [], []
                    for line in option_lines:
                        candidates = [name for name in options if line.startswith('--' + name)]
                        if not candidates:
                            actual_names.append('<unknown>')
                            rendering_errors.append(line)
                            continue
                        name = max(candidates, key=len)
                        actual_names.append(name)
                        # Source pads to width 25, with zero padding for long
                        # IDs: description text can touch the ID directly.
                        prefix = ('--' + name).ljust(25)
                        desc = options[name]['description']
                        if desc is None:
                            if not line.startswith(prefix):
                                rendering_errors.append(line)
                        elif line != prefix + desc.splitlines()[0]:
                            rendering_errors.append(line)
                    record['expected_documentation_options'] = documentation
                    record['actual_documentation_options'] = actual_names
                    record['documentation_rendering_errors'] = rendering_errors
                    passed = passed and actual_names == documentation and not rendering_errors
            except subprocess.TimeoutExpired:
                passed = False
                record = {'kind': kind, 'args': argv, 'failure': 'wall-clock limit exceeded'}
        record['passed'] = passed
        cases.append(record)

    stop_error = 'Error: Invalid option: ' + sentinel
    for name in commands:
        expected = stop_error if name in baseline['commands'] else f'Error: Invalid command: {name}. To print help message: diamond help'
        run('command_registration', [name, '--' + sentinel], expected)
    run('command_leading_double_dash', ['--version', '--' + sentinel], stop_error)

    for name, option in options.items():
        arity = option['arity']
        if arity['kind'] == 'unknown':
            raise ValueError(f'Cannot probe unknown option arity: {name}')
        registered = name in baseline['options']
        group = groups[option['group']]
        disabled = group['disabled_expression'] == 'true'
        allowed = [command for command in baseline['commands']
                   if commands[command]['enum_expression'] in group['command_enum_expressions']]
        command = allowed[0] if allowed else 'blastp'
        values = [] if arity['kind'] == 'bool' else ['1'] * arity['minimum']
        expected = 'Error: Invalid option: ' + name if not registered or disabled else stop_error
        run('long_registration_arity', [command, '--' + name, *values, '--' + sentinel], expected)
        if not registered or disabled:
            continue
        denied = next((c for c in baseline['commands'] if c not in allowed), None)
        if denied:
            run('workflow_permission', [denied, '--' + name, *values, '--' + sentinel],
                'Error: Option is not permitted for this workflow: ' + name)
        short = option['short_id']
        count_error = f"Error: Invalid parameter count for option '{'-' + short + '/' if short else ''}--{name}'"
        if arity['kind'] == 'bool':
            wrong = [['1']]
        elif arity['kind'] == 'scalar':
            wrong = [[], ['1', '2']]
        elif arity['minimum'] > 0:
            wrong = [[]]
        else:
            wrong = []
        for bad in wrong:
            run('arity_rejection', [command, '--' + name, *bad, '--' + sentinel], count_error)
        if arity['kind'] == 'vector':
            run('vector_multiple_values', [command, '--' + name, '1', '2', '--' + sentinel], stop_error)
    for short, name in baseline['short_option_map'].items():
        option = options[name]
        group = groups[option['group']]
        disabled = group['disabled_expression'] == 'true'
        allowed = [c for c in baseline['commands'] if commands[c]['enum_expression'] in group['command_enum_expressions']]
        command = allowed[0] if allowed else 'blastp'
        values = [] if option['arity']['kind'] == 'bool' else ['1'] * option['arity']['minimum']
        expected = 'Error: Invalid option: ' + short if disabled else stop_error
        run('effective_short_registration', [command, '-' + short, *values, '--' + sentinel], expected)
        if values:
            run('effective_short_joined_value', [command, '-' + short + values[0], *values[1:], '--' + sentinel], expected)
    # These four commands print documentation, then necessarily throw on the
    # missing database in Config itself. No command handler can execute.
    for command in ['makedb', 'blastp', 'blastx', 'dbinfo']:
        run('documentation_before_missing_database', [command],
            'Error: Missing parameter: database file (--db/-d)',
            baseline['documentation_options'][command])
    after = sha(binary)
    if before != after:
        raise ValueError('Oracle binary changed during registry probes')
    result = {'oracle_sha256_before': before, 'oracle_sha256_after': after,
              'variant': 'baseline', 'sentinel': sentinel,
              'limits': {'core_bytes': 0, 'cpu_seconds': 2, 'address_space_bytes': 512 << 20,
                         'each_output_file_bytes': cap, 'wall_seconds': 3},
              'source_registry_sha256': artifact['source_sha256'],
              'count': len(cases), 'passed': sum(case['passed'] for case in cases),
              'failed': [case for case in cases if not case['passed']], 'cases': cases,
              'meaning': 'Parser registration, arity, permissions, short-ID resolution and selected documentation only. No alignment/command semantics were executed.'}
    destination.parent.mkdir(parents=True, exist_ok=True)
    destination.write_text(json.dumps(result, indent=2) + '\n')
    print(f"Parser probes: {result['passed']}/{result['count']} pass; evidence: {destination}")
    return not result['failed']


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--repo', type=Path, default=Path(__file__).resolve().parent.parent)
    parser.add_argument('--write', type=Path)
    parser.add_argument('--check', type=Path)
    parser.add_argument('--probe', type=Path, help='Save bounded parser-only oracle probe evidence')
    parser.add_argument('--oracle', type=Path, help='Pinned unchanged baseline executable')
    args = parser.parse_args()
    try:
        artifact = extract(args.repo.resolve())
        serialized = json.dumps(artifact, indent=2, ensure_ascii=False) + '\n'
        if args.write:
            args.write.parent.mkdir(parents=True, exist_ok=True)
            args.write.write_text(serialized)
        if args.check and args.check.read_text() != serialized:
            print('Registry artifact is stale; regenerate with --write.', file=sys.stderr)
            return 1
        print(json.dumps({'commands': {k: v['command_count'] for k, v in artifact['build_variants'].items()}, 'options': {k: {'registered': v['option_count'], 'enabled': v['enabled_option_count']} for k, v in artifact['build_variants'].items()}, 'coverage': artifact['rust_name_coverage'], 'failures': artifact['audit_failures']}, indent=2))
        probes_ok = True
        if args.probe:
            probes_ok = probe_registry(args.repo.resolve(), artifact, args.oracle or args.repo / 'diamond/build/diamond', args.probe)
        return int(bool(artifact['audit_failures']) or not probes_ok)
    except (ValueError, OSError, KeyError, IndexError, SyntaxError, ImportError) as error:
        print(f'Registry audit failed: {error}', file=sys.stderr)
        return 1


if __name__ == '__main__':
    sys.exit(main())
