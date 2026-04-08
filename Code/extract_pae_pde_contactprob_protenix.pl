#!/usr/bin/perl
use strict;
use warnings;
use JSON;

# ✅ Usage: input.json field_name output.csv
die "Usage: $0 protenix.json (token_pair_pae|token_pair_pde|contact_probs) output.csv\n"
    unless @ARGV == 3;

my ($json_file, $field, $outfile) = @ARGV;

# ✅ Validate requested field
my %allowed = map { $_ => 1 } qw(token_pair_pae token_pair_pde contact_probs);
die "❌ Field must be one of: token_pair_pae, token_pair_pde, contact_probs\n"
    unless exists $allowed{$field};

# ✅ Read JSON file
open my $fh, "<", $json_file or die "Can't open $json_file: $!";
local $/;
my $json_text = <$fh>;
close $fh;

# ✅ Parse top-level JSON (should be one model)
my $data = decode_json($json_text);

# ✅ Extract field directly
my $matrix = $data->{$field}
    or die "Field '$field' not found in JSON file.\n";

# ✅ Decode if stored as a string
if (!ref($matrix)) {
    eval {
        $matrix = decode_json($matrix);
    };
    die "Failed to decode '$field' string: $@" if $@;
}

# ✅ Confirm it's a 2D array
die "'$field' is not a 2D matrix\n" unless ref($matrix->[0]) eq 'ARRAY';

# ✅ Write to CSV
open my $out, ">", $outfile or die "Can't write to $outfile: $!";
foreach my $row (@$matrix) {
    print $out join(",", @$row), "\n";
}
close $out;

print "✅ Saved '$field' matrix to: $outfile\n";

