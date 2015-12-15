#!/usr/bin/env python2.7

"""
Quick script to generate Azure Account SAS credentials, to grant people access
to your blob storage containers (per account).

"""

import hmac
import hashlib
import base64
import urllib
import argparse
import sys
import os
import doctest

# What version of SAS are we trying to generate?
SUPPORTED_VERSION="2015-04-05"

# This holds mappings from human-readable names to Account ASA query parameter
# names.
NAME_MAPPINGS = {
    "SignedVersion": "sv",
    "SignedServices": "ss",
    "SignedResourceTypes": "srt",
    # Note that this is singular but can still be a bunch of characters.
    # Docs typo?
    "SignedPermission": "sp", 
    "SignedStart": "st",
    "SignedExpiry": "se",
    "SignedIP": "sip",
    "SignedProtocol": "spr",
    "Signature": "sig"
}

def parse_args(args):
    """
    Takes in the command-line arguments list (args), and returns a nice argparse
    result with fields for all the options.
    
    Borrows heavily from the argparse documentation examples:
    <http://docs.python.org/library/argparse.html>
    """
    
    # Construct the parser (which is stored in parser)
    # Module docstring lives in __doc__
    # See http://python-forum.com/pythonforum/viewtopic.php?f=3&t=36847
    # And a formatter class so our examples in the docstring look good. Isn't it
    # convenient how we already wrapped it to 80 characters?
    # See http://docs.python.org/library/argparse.html#formatter-class
    parser = argparse.ArgumentParser(description=__doc__, 
        formatter_class=argparse.RawDescriptionHelpFormatter)
    
    # General options
    parser.add_argument("--account_name", required=True,
        help="Name of the storage account to operate on")
    parser.add_argument("--account_key", required=True,
        help="Key of the storage account to operate on")
    parser.add_argument("--expiry", required=True,
        help="Expiration date (YYYY-MM-DD)")
        
    # Mode
    parser.add_argument("--mode", choices=["account"], default="account",
        help="Type of SAS togen to generate (account only right now)")
    
    # Account mode options
    parser.add_argument("--services", default="b",
        help="Services to grant access to (of 'b', 'q', 't', and 'f')")
    parser.add_argument("--resource_types", default="co",
        help="Resource types to grant access to (of 's', 'c', and 'o')")
    parser.add_argument("--permissions", default="rl",
        help="Permissions to grant (of 'r', 'w', 'd', 'l', 'a', 'c', 'u', 'p')")

    
    # The command line arguments start with the program name, which we don't
    # want to treat as an argument for argparse. So we remove it.
    args = args[1:]
        
    return parser.parse_args(args)

def construct_string_to_sign(account_name, parameters, mode="account"):
    """
    Compose a UTF-8 encoded SAS string to sign from an account
    name and a dict defining the SAS parameters ("SignedServices",
    "SignedResourceTypes", etc.).
    
    All SAS parameters in the dict must be non-URL-encoded UTF-8-encoded str
    objects, and all required-by-Azure parameters must be filled in.
    
    TODO: Microsoft does not specify how to handle missing optional fields when
    computing the HMAC.
    """
    
    # We'll build up a list of parts.
    parts = []
    
    if mode == "account":
        # We're doing this at the account level.
        
        # From the docs:
        # StringToSign = accountname + "\n" +
        # signedpermissions + "\n" +
        # signedservice + "\n" +
        # signedresourcetype + "\n" +
        # signedstart + "\n" +
        # signedexpiry + "\n" +
        # signedIP + "\n" +
        # signedProtocol + "\n" +
        # signedversion + "\n"
        
        # First in the spec, we have the account name.
        parts.append(account_name)
        
        # Next, the string of permission characters
        parts.append(parameters.get("SignedPermission", ""))
        
        # Next, the service character string
        parts.append(parameters.get("SignedServices", ""))
        
        # Next, the resource type list
        parts.append(parameters.get("SignedResourceTypes", ""))
        
        # Next, the start time
        parts.append(parameters.get("SignedStart", ""))
        
        # Next, the expiration time
        parts.append(parameters.get("SignedExpiry", ""))
        
        # Next, the IP addresses that are allowed
        parts.append(parameters.get("SignedIP", ""))
        
        # Next, the protocols that are allowed (http,https)
        parts.append(parameters.get("SignedProtocol", ""))
        
        # And finally, the version for signature checking.
        parts.append(parameters.get("SignedVersion", ""))
    
    # Add a trailing empty string to get the trailing newline
    parts.append("")
    
    # Join on newlines
    return "\n".join(parts)
    
def sign(account_name, account_key, sas_parameters, mode="account"):
    """
    Sign the given sas parameters for the given account name with the given
    Base64-encoded account key. Returns a Base64-encoded signature string.
    
    """
    
    # Construct the string that needs to be signed, for the given mode
    string_to_sign = construct_string_to_sign(account_name, sas_parameters,
        mode=mode)
    
    print("Signing {}".format(string_to_sign))
    
    # Compute the digest of the SHA256 HMAC with the right key and data
    hmac_digest = hmac.new(base64.b64decode(account_key), string_to_sign,
        hashlib.sha256).digest()
    
    # Return the encoded digest
    return base64.b64encode(hmac_digest)
    
def create_account_sas_parameters(services, resource_types, permissions,
    expiry):
    """
    Create a populated parameters dict for "account" mode from only a services
    string (like "b" for blobs), a resource types string (like "co" for
    container and object access), a permissions string (like "rl" for read and
    list), and an expiration date (like "2016-12-31")
    
    """
    
    # Make the dict. The only special thing we do is set the version to the one
    # we were written against, and set protocols to the specified "default
    # value" in the docs.
    parameters = {
        "SignedVersion": SUPPORTED_VERSION,
        "SignedServices": services,
        "SignedResourceTypes": resource_types,
        "SignedPermission": permissions,
        "SignedExpiry": expiry,
        "SignedProtocol": "https,http"
    }
    
    return parameters
    
def encode_sas_parameters(parameters, signature):
    """
    Given an SAS parameters dict and a Base64-encoded signature, produce the URL
    query string fragment encoding the parameters and signature.
    
    """
    
    # Compose the query parameter key-value mapping
    # Start with all the parameters with their short names
    query_mapping = {NAME_MAPPINGS[k]: v for k, v in parameters.iteritems()}
    
    # Add in the signature
    query_mapping["sig"] = signature
    
    # Encode this dict as a query string and return it.
    return urllib.urlencode(query_mapping)
    
def main(args):
    """
    Parses command line arguments and do the work of the program.
    "args" specifies the program arguments, with args[0] being the executable
    name. The return value should be used as the program's exit code.
    """
    
    if len(args) == 2 and args[1] == "--test":
        # Run the tests
        return doctest.testmod(optionflags=doctest.NORMALIZE_WHITESPACE)
    
    options = parse_args(args) # This holds the nicely-parsed options object
    
    if options.mode == "account":
    
        # Make some account SAS parameters
        sas_parameters = create_account_sas_parameters(options.services,
            options.resource_types, options.permissions, options.expiry)
            
    # Sign them in the correct mode
    signature = sign(options.account_name, options.account_key, sas_parameters,
        mode=options.mode)
        
    # Encode them and print the result
    print("SAS {} parameter string: {}".format(options.mode,
        encode_sas_parameters(sas_parameters, signature)))
    
if __name__ == "__main__" :
    sys.exit(main(sys.argv))
    
    
    
