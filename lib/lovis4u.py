# %%writefile lovis4u_module.py
import lovis4u

def generate_lovis_plot(gff_path, output_pdf_path):
    """
    Generates a lovis4u plot from a GFF file.

    Args:
        gff_path (str): Path to the input GFF file.
        output_pdf_path (str): Path to save the output PDF file.
    """
    # Creating a parameters object and loading config
    parameters = lovis4u.Manager.Parameters()
    parameters.load_config("standard")

    # Example of changing a particular parameter
    parameters.args["output_dir"] = "." # Use current directory for output
    parameters.args["filename"] = output_pdf_path # Set the output filename in parameters

    # To turn off progress messages
    parameters.args["verbose"] = False

    # Creating a loci object and loading gff files
    loci = lovis4u.DataProcessing.Loci(parameters=parameters)

    # Loading folder with gff files
    loci.load_loci_from_extended_gff(gff_path)

    # Set colours  (optional)
    loci.set_feature_colours_based_on_groups()
    loci.set_category_colours()

    # Defining labels to be shown (optional)
    loci.define_labels_to_be_shown()

    # Saving annotation tables (optional)
    loci.save_feature_annotation_table()
    loci.save_locus_annotation_table()

    # Visualisation steps
    # Creating a canvas manager object
    canvas_manager = lovis4u.Manager.CanvasManager(parameters)
    canvas_manager.define_layout(loci)

    # Adding tracks. The only mandatory: loci
    canvas_manager.add_loci_tracks(loci)

    # We can add scale line on the bottom (optional)
    canvas_manager.add_scale_line_track()

    # Category colours (optional)
    canvas_manager.add_categories_colour_legend_track(loci)

    # And homology line track (optional)
    canvas_manager.add_homology_track()

    # Finally, plotting results and saving the pdf file
    canvas_manager.plot(filename=output_pdf_path) # Pass filename explicitly