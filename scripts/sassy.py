#!/usr/bin/env python2.7

"""
Quick script to generate Azure SAS credentials, to grant people access
to your blob storage containers.

This script supports "Account SAS" (--mode account), which lets you grant access
to services of a certain type across all containers in an account, and "Service
SAS" (--mode service), which lets you grant access to a single blob or
container.

Currently, in service mode, only the blob service is supported.

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
    # Fields for account mode
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
    "Signature": "sig",
    # Additional fields for service mode
    "SignedResource": "sr",
    "SignedIdentifier": "si",
    # Parameters for forcing headers on returned blobs
    "Cache-Control": "rscc",
    "Content-Disposition": "rscd",
    "Content-Encoding": "rsce",
    "Content-Language": "rscl",
    "Content-Type": "rsct",
    # Table stuff
    "TableName": "tn",
    "StartPK": "spk",
    "StartRK": "srk",
    "EndPK": "epk",
    "EndRK": "erk"
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
    
    def add_common_options(subparser):
        """
        We want coptions commoin across subparsers, but don't want them before
        the subcommand name.
        """
        
        subparser.add_argument("--account_name", required=True,
            help="Name of the storage account to operate on")
        subparser.add_argument("--account_key", required=True,
            help="Key of the storage account to operate on")
        subparser.add_argument("--expiry", required=True,
            help="Expiration date and time (ISO 8601; YYYY-MM-DD)")
        subparser.add_argument("--permissions", required=True,
            help="Permissions to grant (of 'rwdlacup')")
    
    # Enable subcommands
    subparsers = parser.add_subparsers(
        help="Choose type of SAS token to generate")
        
    # Account mode options
    account_parser = subparsers.add_parser("account",
        help="Generate Account SAS token")
    account_parser.set_defaults(mode="account")
    add_common_options(account_parser)
    account_parser.add_argument("--services", default="b",
        help="Services to grant access to (of 'b', 'q', 't', and 'f')")
    account_parser.add_argument("--resource_types", default="co",
        help="Resource types to grant access to (of 's', 'c', and 'o')")
        
    # Service mode options
    service_parser = subparsers.add_parser("service",
        help="Generate Service SAS token")
    service_parser.set_defaults(mode="service")
    add_common_options(service_parser)
    service_parser.add_argument("--resource_type", default="c",
        help="Type of the resource for the service ('b', 'c', 'f', or 's')")
    service_parser.add_argument("--service", choices=["blob"], default="blob",
        help="Service type (blob only for now, table and queue eventually)")
    service_parser.add_argument("--path", required=True,
        help="path to grant access to (container[/path/to/blob.file])")

    
    # The command line arguments start with the program name, which we don't
    # want to treat as an argument for argparse. So we remove it.
    args = args[1:]
        
    return parser.parse_args(args)

def construct_string_to_sign(account_name, parameters, mode="account",
    service="blob", path=""):
    """
    Compose a UTF-8 encoded SAS string to sign from an account
    name and a dict defining the SAS parameters ("SignedServices",
    "SignedResourceTypes", etc.).
    
    All SAS parameters in the dict must be non-URL-encoded UTF-8-encoded str
    objects, and all required-by-Azure parameters must be filled in.
    
    In "service" mode, the service parameter gives the type of service we're
    granting access to. Currently only "blob" is supported, but individual blobs
    or whole containers of blobs will both work.
    
    In "service" mode, the path parameter must hold the path to the resource to
    which we are granting access, beginning with the container name. So for a
    container it is just the container name, while for a blob it is
    "container/path/to/blob.file".
    
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
        
    elif mode == "service":
        # We're doing this at the service level, for a blob
        
        # Docs say:
        # StringToSign = signedpermissions + "\n" +
        # signedstart + "\n" +
        # signedexpiry + "\n" +
        # canonicalizedresource + "\n" +
        # signedidentifier + "\n" +
        # signedIP + "\n" +
        # signedProtocol + "\n" +
        # signedversion + "\n" +
        # rscc + "\n" +
        # rscd + "\n" +
        # rsce + "\n" +
        # rscl + "\n" +
        # rsct
        
        parts.append(parameters.get("SignedPermission", ""))
        parts.append(parameters.get("SignedStart", ""))
        parts.append(parameters.get("SignedExpiry", ""))
        
        # Add in the canonicalized resource, which we get by tacking the service
        # and account onto the path we're operating on. The docs say there
        # should be no leading slash, but the .NET implementation from Microsoft
        # does include a leading slash.
        # See <https://github.com/Azure/azure-storage-net/blob/40c04ce1dbbc0c569316da0e3fe47b9bff5d3042/Lib/Common/Blob/CloudBlobContainer.Common.cs#L171-L179>
        canonicalized_resource = "/{}/{}/{}".format(service, account_name, path)
        parts.append(canonicalized_resource)
        
        parts.append(parameters.get("SignedIdentifier", ""))
        parts.append(parameters.get("SignedIP", ""))
        parts.append(parameters.get("SignedProtocol", ""))
        parts.append(parameters.get("SignedVersion", ""))
        
        parts.append(parameters.get("Cache-Control", ""))
        parts.append(parameters.get("Content-Disposition", ""))
        parts.append(parameters.get("Content-Encoding", ""))
        parts.append(parameters.get("Content-Language", ""))
        parts.append(parameters.get("Content-Type", ""))
        
        # No trailing newline here
    
    # Join on newlines
    return "\n".join(parts)
    
def sign(account_name, account_key, sas_parameters, mode="account",
    service="blob", path=""):
    """
    Sign the given SAS parameters for the given account name with the given
    Base64-encoded account key. Returns a Base64-encoded signature string.
    
    In "service" mode, "service" specifies the kind of service to grant access
    to, and "path" gives the path to that service/item, beginning with the
    container name.
    
    """
    
    # Construct the string that needs to be signed, for the given mode
    string_to_sign = construct_string_to_sign(account_name, sas_parameters,
        mode=mode, service=service, path=path)
    
    print("Signing:\n--------\n{}\n--------".format(string_to_sign))
    
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
    
def create_blob_service_sas_parameters(resource_type, permissions, expiry):

    """
    Generates a populated SAS parameter dict for "service" mode for a blob
    service, given only the resource type (like "c" for containers or "b" for
    individual blobs), the permissions (from "rwdl" for containers and "racwd"
    for blobs), and the expiry date (like "2016-12-31").
    
    The actual container and blob name (and the fact that this is for blobs)
    aren't specified here (or in the URL). They only factor into the signature.
    
    """
    
    # Make the dict. The only special thing we do is set the version to the one
    # we were written against, and set protocols to the specified "default
    # value" in the docs.
    parameters = {
        "SignedVersion": SUPPORTED_VERSION,
        "SignedResource": resource_type,
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
            
    elif options.mode == "service":

        # We have to make different parameters for each kind of service
        if options.service == "blob":
    
            # Make some service SAS parameters
            sas_parameters = create_blob_service_sas_parameters(
                options.resource_type, options.permissions, options.expiry)
                
        else:
            # No other service types supported yet
            raise RuntimeError("Unsupported SAS service: {}".format(
                options.service))
    
    else:
        # That's all we can do
        raise RuntimeError("Unsupported SAS mode: {}".format(options.mode))
            
    # Sign them in the correct mode. Pass along service and path in case they
    # are defined.
    signature = sign(options.account_name, options.account_key, sas_parameters,
        mode=options.mode, service=getattr(options, "service", None),
        path=getattr(options, "path", None))
        
    # Encode them and print the result. Encoding is mode independent.
    print("SAS {} parameter string: {}".format(options.mode,
        encode_sas_parameters(sas_parameters, signature)))
    
if __name__ == "__main__" :
    sys.exit(main(sys.argv))
    
    
    
