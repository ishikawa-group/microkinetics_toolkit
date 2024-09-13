if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser()
    parser.add_argument("--unique_id")

    args = parser.parse_args()

    unique_id = args.unique_id

    print(f"unique_id = {unique_id}")
