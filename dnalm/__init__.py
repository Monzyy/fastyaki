import os
import sys
from argparse import ArgumentDefaultsHelpFormatter, ArgumentParser


def main():
    parser = ArgumentParser('dnalm', formatter_class=ArgumentDefaultsHelpFormatter)

    subparsers = parser.add_subparsers(
        title='subcommands', description='valid commands',
        help='additional help', dest='command'
    )
    subparsers.required = True

    for module in ('buildlm', ):
        mod = globals()[module]
        p = subparsers.add_parser(module, parents=[mod.argparser()])
        p.set_defaults(func=mod.main)

    args = parser.parse_args()
    args.func(args)


if __name__ == '__main__':
    this_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    sys.path.append(this_dir)
    from dnalm import buildlm
    main()
