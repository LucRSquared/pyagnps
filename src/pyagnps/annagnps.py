def format_mgmt_operation_for_output(df):
    df['Effect_Code_01'] = df['Effect_Code_01'].astype('Int64')
    df['Effect_Code_02'] = df['Effect_Code_02'].astype('Int64')
    df['Effect_Code_03'] = df['Effect_Code_03'].astype('Int64')
    df['Effect_Code_04'] = df['Effect_Code_04'].astype('Int64')
    df['Effect_Code_05'] = df['Effect_Code_05'].astype('Int64')
    return df

def format_mgmt_schedule_for_output(df):
    df['Event_Year'] = df['Event_Year'].astype('Int64')
    df['Event_Month'] = df['Event_Month'].astype('Int64')
    df['Event_Day'] = df['Event_Day'].astype('Int64')

    return df